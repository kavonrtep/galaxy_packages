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
