# Repository guide for coding agents

This is the project instruction set for coding agents working in this repository.
`AGENTS.md` is a symlink to this file, so Claude Code and Codex read the same
text — edit this file, never the symlink. Keeping them as two copies let them
drift: AGENTS.md sat 62 lines behind and was missing the conda channel-order
pitfall, which was live in `~/.planemo.yml` for weeks.

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
  keeps tmp dirs, scratch, and large binary DB files out of the published tarball. When you
  add files that should not ship (scratch data, generated indexes), add them to `exclude:`.
  **Do not exclude the fixtures a `<tests>` block references** — see the rule below.
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
  --conda_channels conda-forge,bioconda,petrnovak \
  --conda_dependency_resolution --conda_auto_install --conda_auto_init \
  --galaxy_root <persistent>/gx             # reuse across runs so Galaxy installs once

planemo serve <tool1>.xml <tool2>.xml --host 127.0.0.1 --port 9090 \
  --galaxy_root <persistent>/gx --galaxy_branch release_25.1 \
  --conda_prefix /home/petr/miniconda3 --conda_channels conda-forge,bioconda,petrnovak \
  --conda_dependency_resolution --conda_auto_install --conda_auto_init
```

Pitfalls hit before:
- **Channel order matters: `conda-forge` must come first.** Galaxy's conda resolver tries
  `--strict-channel-priority` first (`lib/galaxy/tool_util/deps/conda_util.py`), and under strict
  priority a channel list that demotes `conda-forge` below `bioconda` pushes the solver onto
  bioconda's ancient R packages. `tidecluster` then fails to resolve (`r-igraph 2.0.3` ->
  `glpk >=5.0` conflict) on every version, 1.18.0 included. This repo previously documented
  `petrnovak,bioconda,conda-forge`, which is exactly that broken order; `conda-forge,bioconda,petrnovak`
  resolves under both strict and flexible priority. Add `petrnovak` by appending it, never by
  prepending. (Galaxy itself recovers either way — it retries without strict priority — but a
  hand-run `conda create` does not.)
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
- **Run these conda tools with the environment activated, never by absolute path.**
  `/path/to/envs/__dante_tir@0.3.1/bin/dante_tir.py` leaves `cap3`, `mmseqs` and `blastn` off
  `PATH`. The run still **exits 0** and reports zero results, so it is indistinguishable from a
  genuinely negative run; the cause appears only in a per-step file inside the working directory
  (`working_dir/*.cap.err` → `/bin/sh: 1: cap3: not found`). A whole investigation was built on
  unactivated runs and had to be thrown away. Use
  `source $CONDA/bin/activate <env> && dante_tir.py …`. Galaxy activates the env itself, so this
  only bites manual runs — and when a run fails, read the per-step `.err`/`.log` files in the
  working directory before concluding anything from the top-level message, which for these
  pipelines is often just `error in running command`.
- **`dante_tir` needs a channel carrying `r-rbeast`** (every published version depends on it).
  It is absent from conda-forge, bioconda and petrnovak, and comes from `pkgs/r` — part of conda's
  `defaults` — or from the `r` channel. So it resolves under a normal conda config, but a solve
  with `--override-channels -c conda-forge -c bioconda -c petrnovak` fails with `nothing provides
  r-rbeast`, which looks like a broken package and is not. For `planemo test` use
  `--conda_channels conda-forge,bioconda,petrnovak,r` — appended, never prepended, since
  `conda-forge` must stay first.
- **`detect_errors="aggressive"` is unusable for anything that runs mmseqs2.** Aggressive fails a
  job on any stdout/stderr line matching `error:` (case-insensitive, fatal), and mmseqs2 prints
  `there must be an error: N deleted from M that now is empty, but not assigned to a cluster` as
  ordinary clustering chatter — the string is compiled into the binary (seen in 16.747c6). A
  51-minute `tidecluster` run was marked failed on that line alone: exit status 0, all five outputs
  written at full size, archive passing `unzip -t`. Galaxy also marks the output datasets `error`,
  so they cannot be chained onward. `tidecluster/macros.xml` now defines a `stdio` macro that keeps
  the two halves of aggressive that do not misfire — non-zero exit (trustworthy since TideCluster
  1.20.1 made failing steps abort) and the four out-of-memory patterns — and drops the generic
  `error:`/`exception:` matches. The same macro is now duplicated into `dante_ltr`, `dante_tir`
  and `carp` (each is its own Tool Shed repo, so the macro cannot be shared) and expanded by all
  of their tools. `dante_tir/dante_tir.xml` had a worse variant of the same bug: an explicit
  `<stdio>` block with `<regex match="error" ...level="fatal"/>`, which fires on the bare
  substring anywhere in stderr. Note that `<stdio>` rules are *additive* to `detect_errors`
  rather than overriding it (`parse_stdio` in `lib/galaxy/tool_util/parser/xml.py` prepends
  them), so a tool can carry both. The two remaining `aggressive` tools, `dante/summarize_gff.xml`
  and `dante/dante_gff_to_tabular.xml`, only run `summarize_gff.R` over a GFF3 and never reach
  mmseqs2, so they were left alone.

Publishing to the Tool Shed (run from the repo root, pass the tool dir):
- **testtoolshed** (sandbox, failures OK), owner `petrn`:
  `planemo shed_update --shed_target testtoolshed --shed_key $KEY_TEST --owner petrn --force_repository_creation tidecluster/`
  (`--force_repository_creation` both creates the repo and uploads on first push; `shed_create`
  errors out if the repo already exists.)
- **main toolshed** (only after tests pass), owner `petr-novak`: **Petr's manual step — do not
  run it.** He publishes from inside the tool directory with his production key:
  `planemo shed_update --shed_target toolshed --shed_key $KEY --owner petr-novak .`
  Because it happens outside the conversation, **never state what is published from memory** —
  query the Tool Shed API in the same turn you make the claim. See the
  `published-version-check` skill for that and for the GitHub / GHCR / Anaconda equivalents
  (see below for what access is available).
- **What network access you have depends on whether you are sandboxed.** Check, do not assume.
  Outside the sandbox, `gh` is installed and authenticated as `kavonrtep`, so the GitHub API,
  issues and HTTPS pushes all work; ssh egress does not, so `git push` over the
  `git@github.com:` remote fails regardless. Inside the sandbox, assume read-only HTTPS via
  `curl` only. `gh auth status` settles it in one command.
- **Always ask before any outward-facing `gh` action.** Filing issues, commenting and pushing
  happen under Petr's account, so they need his word each time — having the token is not
  standing authorization. Reading (`gh issue view`, `gh api`) is fine unasked.

- **Test fixtures must ship.** A fixture referenced by a `<tests>` block belongs in git and in
  the tarball, or the test can only ever run on the machine that happens to hold the file —
  not from a clone, not from an installed repository. Keep fixtures small enough that this is
  painless (a rule of thumb: under 1 MB; a fixed-seed synthetic or a genome slice is ideal) and
  do not put them in `.gitignore` or in `exclude:`. Issue #5 was exactly this: `carp`'s 199 KB
  `genome_micro.fasta` was in neither, so its test was unrunnable by anyone else.
- `exclude:` is for things no test needs: tmp dirs, generated indexes, and bulk data used only
  by the standalone `test_run*.sh` drivers. Exclude those **by path**, and check the path is
  real — `re_utils` carried entries for a `test_data/` directory that had been renamed to
  `test-data/`, so they matched nothing and 66 MB shipped unnoticed for years.
- Genuinely large reference data (hundreds of MB) should not be in the repository at all; it
  belongs in a conda package, a data manager, or the container image.
- Verify before publishing, since `shed_upload` tars the **working directory**, not git — an
  untracked local file in a tool directory will be uploaded:
  `planemo shed_upload --shed_target testtoolshed --tar_only <tool>/ && tar tzvf shed_upload.tar.gz`

## Generated documentation

`carp/macros.xml` carries four **generated** tokens — `@CLASS_LIST@`,
`@CLASS_LIST_TANDEM@`, `@CLASS_SPECIAL@`, `@VOCABULARY_SOURCE@` — holding the list of
classifications a CARP library header may use. Both carp tools expand them in `<help>`, so the
vocabulary is written once and cannot drift between the two.

They come from CARP's own `classification_vocabulary.yaml`, so **re-run the generator whenever
`@CONTAINER_TAG@` in `carp/macros.xml` moves**:

```
python3 scripts/update_carp_vocabulary_help.py                 # rewrite the tokens
python3 scripts/update_carp_vocabulary_help.py --check         # exit 1 if stale
```

`--check` is the guard against the help describing an older vocabulary than the image the tool
runs; it compares against the pinned tag and does not go looking for newer releases. By default
the vocabulary is read from the matching git tag on GitHub (the image is built from it, and that
is a text download rather than a 4 GB pull); `--from-container` reads the copy inside the image,
and `--vocabulary FILE` uses a local copy when there is no network. Do not edit the blocks between
the `BEGIN GENERATED` / `END GENERATED` markers by hand.

The help also tells the reader where the authoritative list is (the YAML in the image, and the
tool's own HTML report), so a stale block is never the only source.

## Conventions

- Multi-tool repos share one `macros.xml` per directory; add a new tool by dropping its `.xml`
  in the directory and reusing the macros, not by duplicating requirement blocks.
- Keep the `<help>` reStructuredText block and the upstream `homepage_url` accurate — they are
  user-facing on the Tool Shed.
