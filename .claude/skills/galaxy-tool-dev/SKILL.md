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
- **Container tools run the whole `<command>` inside the image** — every binary
  in the command (including collection steps) must exist there. The CARP image
  has no `zip`, so `zip -r` failed after the pipeline; build archives with
  `python3`/`tar` instead.

Interactive inspection (click through the report in a real Galaxy):

```
planemo serve <dir>/*.xml --host 127.0.0.1 --port 9090 \
  --galaxy_root ~/.planemo/galaxy_root \
  --conda_dependency_resolution --conda_auto_install --conda_auto_init
```

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
- `<command detect_errors="aggressive">` (fails on non-zero exit **and** on
  `error:`/`exception:` in stderr) in preference to `<stdio>`.
- Recent `profile=` (~1 year back), version from macro tokens. IUC uses
  `@TOOL_VERSION@+galaxy@VERSION_SUFFIX@`; this repo's existing tools instead use
  `<upstream>.<wrapper>` (e.g. `1.18.0.1`) — match the tool you're editing.
- Escaping: Galaxy params `'$p'` single-quoted; shell vars `\${GALAXY_SLOTS:-1}`;
  Cheetah zero/bool gotchas — guard with `#if str($x).strip()`, compare
  booleans as `str($x) == "true"`.
- `<token>` for Cheetah logic, `<xml>` for element trees — never Cheetah in an
  `<xml>` macro. Add `<citation type="doi">` and an `<xrefs><xref type="bio.tools">`
  where available.
