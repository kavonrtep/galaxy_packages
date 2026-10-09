# Test fixtures

## synthetic_satellite.fasta

Input for the `tidecluster` tests. Generated deterministically by
`make_synthetic_satellite.py`: one clean high-copy 172 bp satellite array of about
86 kb - deliberately above the default 50 kb `min_total_length`, so TAREAN runs -
embedded in random flanks. Small enough to take a full `run_all` through
TideHunter, clustering and TAREAN in a Galaxy test, and clean enough that the
resulting cluster (`TRC_1`), its consensus and the report can all be asserted.

To regenerate:

```
python3 make_synthetic_satellite.py
```

## trc_library.fasta

Input for the `tidecluster_reannotate` ("Annotate Genome") test: the unmutated
172 bp monomer of the array above, with the header `>TRC_1#Satellite/synthetic`.
Written by `make_synthetic_satellite.py` from the same seed as the genome, after
the genome, so adding it left `synthetic_satellite.fasta` byte-identical.

The header is the point of the fixture as much as the sequence: RepeatMasker
takes the string after `#` as the repeat class, the tool copies it into the GFF3
as `Name=`, and the test asserts it. Annotating the genome with it collapses to
one region, `chr_test 24996..111000`, against the array's true 25000..111000.

## comparative_run_a.zip, comparative_run_b.zip

Input for the `tidecluster_comparative` test: two miniature TideCluster run
directories, each holding exactly the six files the comparative analysis reads and
nothing else.

```
tc_consensus_dimer_library.fasta             one consensus dimer per cluster
tc_consensus/consensus_sequences_all.fasta   the per-array dimer pool
tc_clustering.gff3                           array coordinates per cluster
tc_annotation.gff3                           annotated regions
tc_annotation.tsv                            annotation per cluster
tc_tarean/SSRS_summary.csv                   which clusters are SSRs
```

14 clusters per run, 10 of them shared between the two, so the fixture exercises
real cross-sample grouping rather than a degenerate one-cluster-per-family split.

The runs carry the prefix `tc`, not the `tidecluster` this tool suite writes. That
is deliberate: the wrapper reads each run's prefix off its `*_clustering.gff3`
rather than assuming one, and the fixture is what keeps that honest.

Taken from `tests/data/comparative` of TideCluster 1.21.3
(https://github.com/kavonrtep/TideCluster, GPL-3.0), zipped one directory per
archive. Upstream derived it from a real *Solanum lycopersicum* run by keeping at
most 12 sequences per cluster, dropping sequences over 2500 bp and dropping
clusters left with fewer than 2 sequences. 25 kB per archive.

To regenerate from a newer TideCluster release, replace `1.21.3` below:

```
base=https://raw.githubusercontent.com/kavonrtep/TideCluster/1.21.3/tests/data/comparative
for s in a b; do
    mkdir -p sample_$s/tc_consensus sample_$s/tc_tarean
    for f in tc_annotation.gff3 tc_annotation.tsv tc_clustering.gff3 \
             tc_consensus_dimer_library.fasta \
             tc_consensus/consensus_sequences_all.fasta \
             tc_tarean/SSRS_summary.csv; do
        curl -sL "$base/sample_$s/$f" -o "sample_$s/$f"
    done
    ( cd sample_$s && zip -qrX ../comparative_run_$s.zip . )
done
```

## annotation_run.zip, annotation_library.fasta

Input for the `tidecluster_annotation` test. Built by
`make_annotation_fixture.py` out of `comparative_run_a.zip`, so this directory
holds no third copy of a real run:

```
python3 make_annotation_fixture.py
```

`annotation_run.zip` (2.9 kB) is a run directory under the prefix the wrapper
passes (`-pr tidecluster`):

```
tidecluster_clustering.gff3                  6 regions, clusters 26, 28 and 36
tidecluster_consensus/TRC_26_dimers.fasta    4 sequences
tidecluster_consensus/TRC_28_dimers.fasta    3 sequences
tidecluster_consensus/TRC_36_dimers.fasta    2 sequences
```

The `<prefix>_consensus/TRC*dimers.fasta` layout is what
`TideCluster.py annotation` globs for, and it is why the comparative fixtures
cannot be reused directly: their clustering GFF3 carries no
`consensus_sequence` attributes either, so neither of the two paths through the
annotation step is available from them.

`annotation_library.fasta` holds two references, each the first dimer of one of
those clusters relabelled in RepeatMasker format - `testsat_A#Satellite/testfam`
and `testsat_B#rDNA/45S`. A reference taken from a cluster's own pool is a
near-exact match for its siblings, which makes the expected annotation
deterministic. TRC_36 is left without a reference on purpose: the test asserts
it comes back as `annotation=NA`, and the two invented class names come only
from the FASTA headers, which is what pins the documented behaviour that the
string after `#` is what gets reported.

The references are per-array dimers rather than TAREAN consensus sequences
because the dimer library in `comparative_run_a.zip` covers only TRC_1, 4, 6, 8
and 11, all of them larger clusters.
