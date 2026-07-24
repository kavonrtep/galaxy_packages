---
name: galaxy-tool-dev
description: >-
  Use when developing, updating, testing, or publishing the Galaxy tool
  wrappers in this repo (galaxy_packages) — editing tool XML / macros.xml /
  .shed.yml, bumping a tool to a new upstream version, running planemo
  lint/test/serve, or pushing to the Tool Shed or testtoolshed. Triggers on
  "planemo", "toolshed", "testtoolshed", "shed_update", "shed_create",
  ".shed.yml", "galaxy tool", "tool wrapper", or a tool version bump.
---

# Galaxy tool development, testing and publishing

Each top-level directory is one Tool Shed repository (see the repo `CLAUDE.md`
for architecture). Most in-house tools (dante*, tidecluster, …) are thin
wrappers around a conda package on the **`petrnovak`** Anaconda channel.

`~/.planemo.yml` already pins `galaxy_branch: release_25.1`, `conda_prefix`,
`conda_ensure_channels: petrnovak,bioconda,conda-forge`, and the testtoolshed
key — so those flags are not repeated below. planemo is installed persistently
via pipx (`~/.local/bin/planemo`).

## 1. Version bump (thin wrappers)

In the tool's `macros.xml`: set `@REQUIREMENT_VERSION@` to the **exact** conda
package version (confirm it exists: `conda search -c petrnovak <pkg>`), and
`@TOOL_VERSION@` to that plus a wrapper-revision digit (e.g. package `1.18.0`
→ tool `1.18.0.1`). Tools without a `macros.xml` carry both inline on `<tool>`
and `<requirement>` — keep them in sync.

## 2. Check what changed upstream

Don't assume the CLI/outputs are stable across a multi-version jump. Clone the
upstream repo at the target tag and diff the argparse and the produced output
filenames against what the wrapper's `<command>` passes and copies. HTML
reports in particular get restructured (e.g. TideCluster 1.9 moved to a v2
`<prefix>_index.html` + `<prefix>_report/` layout with the old pages under
`<prefix>_report_legacy/`); the wrapper must copy the current report tree into
`extra_files_path` or the report renders broken.

## 3. Verify

Fast smoke of the actual behaviour first (build a scratch conda env with the
target version, run the tool on a tiny input, confirm real output filenames),
then planemo:

```
planemo lint <dir>/<tool>.xml
planemo test <dir>/<tool>.xml \
  --conda_dependency_resolution --conda_auto_install --conda_auto_init \
  --galaxy_root ~/.planemo/galaxy_root
```

Add a `<tests>` block if the tool has none. Pipeline output is not
byte-deterministic — assert on stable content markers (`has_text`), not exact
files. Keep large test genomes in `test-data/`, gitignored, and listed in
`.shed.yml` `exclude:`.

Interactive inspection (click through the report in a real Galaxy):

```
planemo serve <dir>/*.xml --host 127.0.0.1 --port 9090 \
  --galaxy_root ~/.planemo/galaxy_root \
  --conda_dependency_resolution --conda_auto_install --conda_auto_init
```

## 4. Publish

testtoolshed first (sandbox, failures OK, owner `petrn`); main toolshed only
after tests pass (owner `petr-novak`, production key passed at run time):

```
planemo shed_update --shed_target testtoolshed --owner petrn \
  --force_repository_creation <dir>/
planemo shed_update --shed_target toolshed --owner petr-novak \
  --shed_key $KEY <dir>/
```

`--force_repository_creation` both creates the repo and uploads on first push;
`shed_create` errors if the repo already exists.

## Pitfalls (all hit in practice)

- Default Galaxy branch is unstable `master` — the pinned `release_25.1` in
  `~/.planemo.yml` avoids it.
- A stale `~/.planemo/gx_venv*` with a dangling python symlink breaks the
  framework install (`Broken symlink … python3`) — `rm -rf` it and rerun.
- SQLite "database is locked" kills the job handler when a slow `conda create`
  runs inside the first job — pre-build the conda env, or just rerun (the
  resolved `__<tool>@<version>` env is cached).
- `planemo serve` leaves a detached `gunicorn` master bound to the port after
  the wrapper is killed — stop it by the PID listening on the port.
- `planemo shed_lint` can hang for minutes — timeout or skip it.
