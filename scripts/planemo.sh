#!/bin/bash
# planemo.sh — one invocation for lint / test / serve on this repo's tool dirs.
#
#   scripts/planemo.sh lint  tidecluster/            [extra planemo args...]
#   scripts/planemo.sh test  tidecluster/            [extra planemo args...]
#   scripts/planemo.sh test  dante_tir/dante_tir.xml
#   scripts/planemo.sh serve tidecluster/            [extra planemo args...]
#
# What it encodes, all of it learned the hard way:
#
#   * Channel order: conda-forge first, `r` appended. Under
#     --strict-channel-priority a list that demotes conda-forge below bioconda
#     pushes the solver onto bioconda's old R stack and tidecluster stops
#     resolving; `r` carries r-rbeast, which every dante_tir release needs.
#   * Separate galaxy_root per mode. test and serve against the same root fight
#     over one SQLite database ("database is locked"), so serve gets gx/ and test
#     gets gx_test/ and the two can run at once.
#   * Test reports are written outside the repository. planemo drops
#     tool_test_output.{html,json} into the working directory, and shed_upload
#     tars the working directory rather than what git tracks, so a report left in
#     a tool directory gets published.
#   * TMPDIR on the big disk: / is routinely near full, /mnt/ssd is not.
#   * serve refuses to start when something already holds the port, instead of
#     failing obscurely. planemo serve leaks a detached gunicorn when the wrapper
#     is killed, and one was found still holding the port — and the shared
#     galaxy_root — six weeks later. The trap below stops that recurring.
set -euo pipefail

MODE=${1:-}
TARGET=${2:-}
if [[ -z $MODE || -z $TARGET ]]; then
    sed -n '2,8p' "${BASH_SOURCE[0]}" >&2
    exit 2
fi
shift 2

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PLANEMO=${PLANEMO:-$HOME/.local/bin/planemo}
CONDA_PREFIX_DIR=${CONDA_PREFIX_DIR:-/home/petr/miniconda3}
CHANNELS=${CHANNELS:-conda-forge,bioconda,petrnovak,r}
BRANCH=${GALAXY_BRANCH:-release_25.1}
WORK=${PLANEMO_WORK:-/mnt/ssd/carp_planemo}
GX_SERVE=$WORK/gx
GX_TEST=$WORK/gx_test
PORT=${PORT:-9090}

export TMPDIR=${TMPDIR:-$WORK/tmp_test}
mkdir -p "$TMPDIR"
cd "$REPO"

CONDA_ARGS=(
    --conda_prefix "$CONDA_PREFIX_DIR"
    --conda_channels "$CHANNELS"
    --conda_dependency_resolution --conda_auto_install --conda_auto_init
)

# PIDs listening on $PORT, or nothing.
port_pids() { ss -ltnpH "sport = :$PORT" 2>/dev/null | grep -oE 'pid=[0-9]+' | cut -d= -f2 | sort -u; }

case $MODE in
lint)
    exec "$PLANEMO" lint "$TARGET" "$@"
    ;;

test)
    label=$(basename "${TARGET%/}" .xml)
    out=$WORK/test_${label}_$(date +%Y%m%d_%H%M%S)
    echo "planemo test $TARGET  (report: $out.html)" >&2
    set +e
    "$PLANEMO" test "$TARGET" \
        --galaxy_root "$GX_TEST" --galaxy_branch "$BRANCH" \
        "${CONDA_ARGS[@]}" \
        --test_output "$out.html" --test_output_json "$out.json" "$@"
    rc=$?
    set -e
    # Belt and braces: planemo still drops these in cwd on some paths, and an
    # untracked file inside a tool directory would be published by shed_upload.
    find "$REPO" -maxdepth 2 -name 'tool_test_output.html' -o -maxdepth 2 -name 'tool_test_output.json' \
        | xargs -r rm -f
    echo "exit $rc  (report: $out.html)" >&2
    exit $rc
    ;;

serve)
    pids=$(port_pids || true)
    if [[ -n ${pids:-} ]]; then
        echo "port $PORT is already in use by PID(s): $pids" >&2
        ps -o pid,lstart,args -p $pids 2>/dev/null | cut -c1-160 >&2
        echo >&2
        echo "A leaked planemo serve looks exactly like this. If it is one:" >&2
        echo "    kill $pids" >&2
        echo "Then rerun, or set PORT=<other> to leave it alone." >&2
        exit 1
    fi
    # Reap the detached gunicorn planemo leaves behind when this script stops.
    cleanup() {
        local p
        p=$(port_pids || true)
        [[ -n ${p:-} ]] && { echo "stopping leftover Galaxy on port $PORT (PID $p)" >&2; kill $p 2>/dev/null || true; }
    }
    trap cleanup EXIT INT TERM
    echo "planemo serve $TARGET on http://127.0.0.1:$PORT" >&2
    "$PLANEMO" serve "$TARGET" --host 127.0.0.1 --port "$PORT" \
        --galaxy_root "$GX_SERVE" --galaxy_branch "$BRANCH" \
        "${CONDA_ARGS[@]}" "$@"
    ;;

*)
    echo "unknown mode: $MODE (expected lint, test or serve)" >&2
    exit 2
    ;;
esac
