# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A collection of **Galaxy Tool Shed repositories** for repetitive-DNA / transposable-element
annotation (the RepeatExplorer2 tool suite and related utilities). This repo contains **tool
wrappers**, not the analysis code itself. Each top-level directory is an independent Tool Shed
repository, published to the Galaxy Tool Shed under owner `petr-novak` (a few under
`repeatexplorer`).

There is no build step and no repo-wide test runner. Work is editing Galaxy tool XML (and, for
a few tools, the bundled Python/R scripts).

## Layout

Each top-level directory = one Tool Shed repository, identified by its `.shed.yml`:

- `dante/`, `dante_ltr/`, `dante_tir/`, `tidecluster/`, `repeatexplorer2/` — **thin wrappers**.
  The actual program is a separate conda package (bioconda / conda-forge / r / the `petrnovak`
  conda channel) and a separate upstream repo (e.g. `github.com/kavonrtep/dante_tir`). The XML
  only declares a `<requirement type="package">` and builds the command line.
- `re_utils/`, `various_galaxy_tools/`, `short_read_simulator/`, `krona/` — **self-contained
  tools**. The Python/R/shell script lives next to its `.xml` in the same directory and is
  shipped inside the repo.
- `repeat_annotation_pipeline/` — a **git submodule** (`github.com/kavonrtep/repeat_annotation_pipeline`).
- `repex_tarean_old/` — legacy, generally leave alone.

## Anatomy of a tool

- `.shed.yml` — Tool Shed metadata: `name`, `owner`, `categories`, and an `exclude:` list that
  keeps tmp dirs, test data, and large binary DB files out of the published tarball. When you
  add files that should not ship (scratch data, generated indexes), add them to `exclude:`.
- `<tool_name>.xml` — the wrapper: `<command><![CDATA[ ... ]]></command>` chains CLI calls with
  `&&`, using Cheetah templating (`#if`, `$param`, `${output}`). Use `\${GALAXY_SLOTS:-1}` for
  CPU count. Self-contained tools invoke their sibling script by name (it must be executable and
  on `$PATH` at install time).
- `macros.xml` (where present) — defines version tokens and the `requirements` macro; imported
  via `<macros><import>macros.xml</import></macros>` and `<expand macro="requirements"/>`.
- `test-data/` or `test_data/` — inputs and expected outputs referenced by `<tests>` blocks.

## Versioning convention (important, easy to get wrong)

The Galaxy **tool version** is the underlying program version plus one extra suffix digit for
wrapper revisions. The **requirement version** is the exact conda package version.

Example (`dante_tir/macros.xml`): program `dante_tir` is at `0.2.0`, so
`@REQUIREMENT_VERSION@ = 0.2.0` (must match the conda package) and `@TOOL_VERSION@ = 0.2.0.1`
(the `.1` is the wrapper revision — bump it when you change the XML without changing the program).

For thin wrappers, bumping a tool to a new upstream release means: update the conda package
version in `@REQUIREMENT_VERSION@` (or the inline `<requirement>` version), then reset/set
`@TOOL_VERSION@` accordingly. Tools without a `macros.xml` (e.g. `dante/dante.xml`) carry the
version and requirement inline on the `<tool>` and `<requirement>` elements — keep those in sync.

## Working with tools

- `re_utils` also has standalone shell drivers (`test_run1.sh`, `test_run2.sh`) that exercise the
  bundled scripts directly, mirroring the Galaxy `<tests>`.
- Conda dependencies: thin-wrapper programs are installed by Galaxy from the declared
  `<requirement>` at package version — do not `pip install`/`mamba install` locally to "make a
  tool work"; fix the requirement declaration instead. The upstream conda packages for the
  in-house tools (dante, dante_ltr, dante_tir, tidecluster, …) live on the **`petrnovak`**
  Anaconda channel, so planemo dependency resolution must include it.

## Testing and publishing with planemo

planemo is the standard tooling but is **not** pinned in this repo — install it once in a
persistent location (`pipx install planemo`, or a dedicated venv; do not rely on a scratch venv
that vanishes between sessions). Conda base is `/home/petr/miniconda3`.

Lint is fast; a full `test` runs the real pipeline through a local Galaxy and takes minutes.
The verified invocation (learned the hard way — see pitfalls):

```
planemo lint <tool>.xml

planemo test <tool>.xml \
  --galaxy_branch release_25.1 \            # pin a release; the default pulls master (unstable)
  --conda_prefix /home/petr/miniconda3 \
  --conda_channels petrnovak,bioconda,conda-forge \
  --conda_dependency_resolution --conda_auto_install --conda_auto_init \
  --galaxy_root <persistent>/gx             # reuse across runs so Galaxy installs once

planemo serve <tool1>.xml <tool2>.xml --host 127.0.0.1 --port 9090 \
  --galaxy_root <persistent>/gx --galaxy_branch release_25.1 \
  --conda_prefix /home/petr/miniconda3 --conda_channels petrnovak,bioconda,conda-forge \
  --conda_dependency_resolution --conda_auto_install --conda_auto_init
```

Pitfalls hit before:
- **Galaxy defaults to `master`** and can fail to build — always pass `--galaxy_branch` (e.g.
  `release_25.1`).
- **Stale `~/.planemo/gx_venv*`** with a dangling python symlink breaks the framework install
  (`Broken symlink … python3`). Fix: `rm -rf` the offending venv and let planemo recreate it.
- **SQLite "database is locked"** kills the job handler when a slow `conda create` runs inside the
  job window (first test of a new tool version). Pre-build the conda env once, or just re-run —
  the resolved `__<tool>@<version>` env is cached and the rerun starts the job immediately.
- **`planemo serve` leaves Galaxy running** after the wrapper is killed: the detached `gunicorn`
  master (bound to the port) survives. Kill it by the PID listening on the port, not by pattern.
- `planemo shed_lint` can hang for minutes — use a timeout or skip it.

Publishing to the Tool Shed (run from the repo root, pass the tool dir):
- **testtoolshed** (sandbox, failures OK), owner `petrn`:
  `planemo shed_update --shed_target testtoolshed --shed_key $KEY_TEST --owner petrn --force_repository_creation tidecluster/`
  (`--force_repository_creation` both creates the repo and uploads on first push; `shed_create`
  errors out if the repo already exists.)
- **main toolshed** (only after tests pass), owner `petr-novak`, `$KEY` supplied at runtime:
  `planemo shed_update --shed_target toolshed --shed_key $KEY --owner petr-novak tidecluster/`
- `.shed.yml` `exclude:` keeps `test-data/` and scratch out of the uploaded tarball; large test
  genomes are also gitignored (kept local).

## Conventions

- Multi-tool repos share one `macros.xml` per directory; add a new tool by dropping its `.xml`
  in the directory and reusing the macros, not by duplicating requirement blocks.
- Keep the `<help>` reStructuredText block and the upstream `homepage_url` accurate — they are
  user-facing on the Tool Shed.
