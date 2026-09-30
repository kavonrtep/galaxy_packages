# Test fixtures

Inputs for the `<tests>` blocks of `dante_tir.xml` and `dante_tir_summary.xml`.
Before these existed, both test blocks referenced eight files that were never in
the repository, so neither tool had a runnable test.

```
test_genome.fasta     1.97 MB   80 concatenated regions of tiny_pea
test_dante.gff3       0.48 MB   640 DANTE records, 174 MuDR/Mutator Subclass_1 copies
test_tir_final.gff3   2.6 kB    DANTE_TIR output for the two above, 4 TIR records
```

`test_tir_final.gff3` is the `DANTE_TIR_final.gff3` that `dante_tir.py 0.3.1`
produces from the other two, so the summary tool's input is consistent with the
genome it is given.

## Why the fixture is 2.4 MB and not smaller

DANTE_TIR assembles the termini of many copies of one TIR superfamily with CAP3,
so it detects nothing until the input holds enough copies. Measured with 0.3.1 on
MuDR/Mutator regions of tiny_pea:

| Subclass_1 copies | fixture size | TIRs detected |
|---|---|---|
| 129 | 1.68 MB | 0 |
| **174** | **2.40 MB** | **4** |
| 204 | 2.93 MB | 5 |
| 244 | 6.80 MB | 7 |

174 copies is the smallest tested input that detects anything, which is why the
fixture sits above the 1 MB rule of thumb in `CLAUDE.md`. The alternatives are
worse: upstream's own `smoke` dataset (2.37 MB) is deliberately below the
threshold and detects nothing, and their `short` dataset that does detect is
20 MB of sequence.

A contiguous slice cannot reach 174 copies at a shippable size — the densest
MuDR/Mutator window in tiny_pea holds just 24 of them in 346 kb. So
`make_test_fixture.py` takes the neighbourhood of each domain (±7 kb) and
concatenates them as separate records, which is roughly 8x more copies per byte.
Regions are chosen densest-first, so `REGIONS` controls size and copy count
together.

The tests assert well below 4 TIRs, so a change in detection sensitivity does not
fail them spuriously — the same reasoning behind upstream's `>= 5` assertion on a
fixture that produced 14.

## Regenerating

```
python3 make_test_fixture.py [SOURCE_DIR] [OUT_DIR]
```

The source dataset is in neither this repository nor upstream: 81 MB FASTA +
22 MB GFF3 at `/mnt/ssd/dante_tir/test-data/` on the development machine. The
script therefore records how the fixtures were made rather than offering to
remake them from a clone. To refresh `test_tir_final.gff3` as well, run
`dante_tir.py` on the regenerated pair and copy its `DANTE_TIR_final.gff3` here.

## Two traps worth knowing

**Activate the conda environment.** Running `dante_tir.py` by absolute path
leaves `cap3`, `mmseqs` and `blastn` off `PATH`. The run still **exits 0** and
reports zero TIRs, so it looks like a fixture that is simply too small; the real
cause appears only in `<output_dir>/working_dir/*.cap.err`
(`/bin/sh: 1: cap3: not found`). Every conclusion drawn from an unactivated run
is worthless. Galaxy activates the environment itself, so this only bites manual
runs.

**`r-rbeast` needs a channel that carries it.** Every published `dante_tir`
depends on it. It is not in conda-forge, bioconda or petrnovak; it comes from
`pkgs/r`, part of conda's `defaults`, or from the `r` channel. It therefore
resolves under a normal conda configuration, but a solve with
`--override-channels -c conda-forge -c bioconda -c petrnovak` fails with
`nothing provides r-rbeast`. For `planemo test`, append `r` to the channel list —
append, never prepend, since `conda-forge` must stay first.

## Version-dependent behaviour the tests rely on

A run that detects no TIRs is handled differently across releases: `0.2.x` exits
0 and writes **no output files at all**, which makes this wrapper's unconditional
`mv output_dir/DANTE_TIR_final.gff3 …` fail on a genome with no detectable TIRs.
`0.3.1` writes all three files, empty, with an explicit message. The wrapper
requires `0.3.1`, so that failure mode is gone; reverting the requirement to
`0.2.x` would bring it back and need a guard in the command.
