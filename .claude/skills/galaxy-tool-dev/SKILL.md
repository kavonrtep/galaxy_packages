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

Each top-level directory is one Tool Shed repository (see the repo instruction
file, `CLAUDE.md` / `AGENTS.md`, for architecture). Most in-house tools (dante*, tidecluster, …) are thin
wrappers around a conda package on the **`petrnovak`** Anaconda channel.

`~/.planemo.yml` pins `galaxy_branch: release_25.1`, `conda_prefix`,
`conda_ensure_channels: conda-forge,bioconda,petrnovak,r`, and the testtoolshed
key. planemo is installed persistently via pipx (`~/.local/bin/planemo`).

**Prefer `scripts/planemo.sh {lint|test|serve} <target>`** over calling planemo
directly. It encodes the channel order, a separate `galaxy_root` per mode, test
reports written outside the repository, `TMPDIR` on the big disk, and a port
check plus cleanup for `serve` — each of which has gone wrong at least once. The
flags below are what it passes, for when you need to deviate.

Channel order matters: `conda-forge` must come **first**, and new channels are
appended, never prepended. Under `--strict-channel-priority` a list that demotes
conda-forge below bioconda pushes the solver onto bioconda's old R stack and
`tidecluster` stops resolving at any version. `r` is in the list because it
carries `r-rbeast`, which every `dante_tir` release depends on and which
conda-forge, bioconda and petrnovak do not have.

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
scripts/planemo.sh lint <dir>/<tool>.xml
scripts/planemo.sh test <dir>/<tool>.xml
```

Use a **different `galaxy_root` for `test` than for `serve`** (the script uses
`gx_test/` and `gx/`). Both against one root means two Galaxy instances on one
SQLite database, which is the "database is locked" failure; with separate roots a
test can run while a serve session stays up.

Add a `<tests>` block if the tool has none. Key rules (learned the hard way):

- **Test at the tool's DEFAULT parameters.** A test that overrides parameters
  can pass while the default path a real user hits is broken — e.g. a synthetic
  fixture whose satellite array was below the default `min_total_length` made
  TAREAN skip the cluster, so the library was empty at defaults while the test
  (run with `-M 5000`) still passed. Size fixtures for the defaults.
- **Assert output *completeness*, not just presence.** A produced-but-empty
  dataset still counts toward `expect_num_outputs`, so a broken output passes
  unless you assert its content: `has_text`, `has_size min="1"`, and for the
  results archive `has_archive_member path="..."` to prove expected files are
  actually inside it.
- Pipeline output is not byte-deterministic — assert on stable content markers,
  not golden files. Use a deterministic fixture (fixed-seed generator) where
  possible.
- **Ship the test data**: commit it under `test-data/` and do NOT `.shed.yml`
  `exclude:` it, so the test runs from a clone and on the Tool Shed. Keep it
  small (< 1 MB; a fixed-seed synthetic is ideal).
- **A detection tool has a threshold, so measure it before sizing a fixture.**
  Tools that assemble evidence across many copies of a repeat detect *nothing*
  below some copy number, and a fixture under it produces a clean, passing,
  completely uninformative run. Build fixtures at two or three sizes, record what
  each detects, and put the table in `test-data/README.md` — `dante_tir` needed
  174 Subclass_1 copies (2.4 MB) where 129 (1.68 MB) found nothing, which is why
  that one fixture is over the 1 MB guideline. Then assert *well below* the count
  you measured, so a shift in detection sensitivity does not fail the test.
  Concatenating only the neighbourhoods of interest as separate records gets far
  more copies per byte than one contiguous slice.
- **Container tools run the whole `<command>` inside the image** — every binary
  in the command (including collection steps) must exist there. The CARP image
  has no `zip`, so `zip -r` failed after the pipeline; build archives with
  `python3`/`tar` instead.

Interactive inspection (click through the report in a real Galaxy):

```
scripts/planemo.sh serve <dir>/
```

It refuses to start when something already holds the port and names the PID,
rather than failing obscurely, and kills the Galaxy it started when you stop it.
Set `PORT=<other>` to run a second session alongside the first.

## 4. Publish

Push to **testtoolshed** (sandbox, failures OK, owner `petrn`; key is in
`~/.planemo.yml`). `--force_repository_creation` creates the repo and uploads
on first push (`shed_create` errors if the repo already exists):

```
planemo shed_update --shed_target testtoolshed --owner petrn \
  --force_repository_creation <dir>/
```

The **main toolshed push is Petr's job — do NOT run it.** Only prepare and
report the command for him to run manually from the tool directory:

```
planemo shed_update --shed_target toolshed --shed_key $KEY --owner petr-novak .
```

## Running a packaged tool by hand

Before blaming a wrapper, reproduce outside Galaxy — but **activate the conda
environment**; do not call the entry point by absolute path:

```
source /home/petr/miniconda3/bin/activate '/home/petr/miniconda3/envs/__dante_tir@0.3.1'
dante_tir.py -g … -f … -o out -c 8
```

Running `<env>/bin/dante_tir.py` directly leaves `cap3`, `mmseqs` and `blastn`
off `PATH`. The run then **exits 0 and reports zero results** — indistinguishable
from a genuinely negative run — with the cause only in
`out/working_dir/*.cap.err` (`/bin/sh: 1: cap3: not found`). A whole
investigation was built on unactivated runs and had to be discarded, including a
wrong report that the conda package was broken. Galaxy activates the environment
itself, so this only bites manual runs.

**When a run fails, read the per-step logs before concluding anything.** These
pipelines deliberately redirect each step into its own `.log`/`.err`, so the
top-level message is often just `error in running command` (R's `system()`) or
snakemake naming a rule. The cause is in `working_dir/*.err`, `out/log/`, or the
rule's own log — never in the summary line. Equally, do not conclude a package is
unavailable from a failed `conda create`: `curl` exit 6 / HTTP 000 and
`repo.anaconda.com` errors are transient DNS, and "fixing" one by narrowing the
channel set manufactures a convincing but false `nothing provides <dep>`.

## Pitfalls (all hit in practice)

- **`--` is illegal inside an XML comment, and tool XML is full of CLI flags.**
  A comment mentioning `--long`, `--cleanenv` or `--min_coverage` makes the file
  unparseable, and `planemo lint` reports it as
  `galaxy.util ERROR: Error parsing file <path>` — naming the file but not the
  cause, which reads like a lint bug rather than a typo. Reword the flag out of
  the comment (`long mode`, `min_coverage of 3`). Hit three times in one session.
  Sweep the whole repo for it with:

  ```
  python3 -c "
  import glob, re, io
  for f in glob.glob('*/[a-z_]*.xml'):
      for m in re.finditer(r'<!--(.*?)-->', io.open(f).read(), re.S):
          if '--' in m.group(1): print(f, m.group(1).strip()[:60])"
  ```

- Default Galaxy branch is unstable `master` — the pinned `release_25.1` in
  `~/.planemo.yml` avoids it.
- A stale `~/.planemo/gx_venv*` with a dangling python symlink breaks the
  framework install (`Broken symlink … python3`) — `rm -rf` it and rerun.
- SQLite "database is locked" kills the job handler when a slow `conda create`
  runs inside the first job — pre-build the conda env, or just rerun (the
  resolved `__<tool>@<version>` env is cached).
- `planemo serve` leaves a detached `gunicorn` master bound to the port after
  the wrapper is killed — stop it by the PID listening on the port, not by
  pattern. **`pkill -f <pattern>` matches your own shell** whenever the pattern
  appears in the command you are running, so `pkill -f "planemo test"` typed
  inside a command containing that string kills the command issuing it, part way
  through. That happened three times in one session, once aborting a cleanup
  during a disk emergency. The safe form collects PIDs first and excludes self:
  `P=$(pgrep -f planemo | grep -vw $$ | grep -vw $PPID); kill $P`.
  **Check the port before serving**: one leak was found still holding
  both the port and the shared `galaxy_root` six weeks later, and a second
  instance on that root would have deadlocked on its database.
  `scripts/planemo.sh serve` checks, and reaps its own on exit.
- A failed job's outputs are marked `error` in Galaxy even when the files are
  complete, so they cannot be used as input downstream. Nothing short of fixing
  the detection and re-running makes them usable again.
- `planemo shed_lint` can hang for minutes — timeout or skip it.

## References & IUC best practices

Consult these for anything non-obvious; they are good agentic resources:

- galaxy-skills tool-dev SKILL: <https://github.com/galaxyproject/galaxy-skills/blob/main/tool-dev/SKILL.md>
- IUC standards / best practices: <https://galaxy-iuc-standards.readthedocs.io/en/latest/best_practices.html>

Key conventions from them, worth applying to new/edited tools here:

- **Element order** (`planemo lint` `XMLOrder`): `description` → `macros` →
  `xrefs` → `requirements` → `stdio`/`version_command` → `command` → `inputs`
  → `outputs` → `tests` → `help` → `citations`. Several older tools here put
  `description` after `macros` — fix when touching them. Run `planemo format`
  before committing.
- IUC suggests `<command detect_errors="aggressive">` over `<stdio>`. **Do not use
  it for anything in this repo whose program reaches mmseqs2** — which is every
  in-house tool. Aggressive fails a job on any stream line matching `error:`, and
  mmseqs2 prints `there must be an error: N deleted from M that now is empty, but
  not assigned to a cluster` as ordinary clustering output. That alone failed a
  completed 51-minute TideCluster run: exit status 0, every output written in
  full, and Galaxy marks the output datasets `error` too, so they cannot be
  chained onward. Expand each repo's `stdio` macro instead (keeps non-zero exit
  and the out-of-memory patterns, drops the text matching). Note also that
  `<stdio>` rules are *additive* to `detect_errors` rather than overriding it, so
  a tool can carry both and the laxer `detect_errors` will not protect it.
- Recent `profile=` (~1 year back), version from macro tokens. IUC uses
  `@TOOL_VERSION@+galaxy@VERSION_SUFFIX@`; this repo's existing tools instead use
  `<upstream>.<wrapper>` (e.g. `1.18.0.1`) — match the tool you're editing.
- Escaping: Galaxy params `'$p'` single-quoted; shell vars `\${GALAXY_SLOTS:-1}`;
  Cheetah zero/bool gotchas — guard with `#if str($x).strip()`, compare
  booleans as `str($x) == "true"`.
- `<token>` for Cheetah logic, `<xml>` for element trees — never Cheetah in an
  `<xml>` macro. Add `<citation type="doi">` and an `<xrefs><xref type="bio.tools">`
  where available.
