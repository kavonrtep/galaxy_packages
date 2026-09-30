#!/usr/bin/env python3
"""Build the dante_tir test fixtures from tiny_pea.

DANTE_TIR assembles the termini of many copies of one TIR superfamily with CAP3,
so it detects nothing until the input holds enough copies. Upstream's smoke
dataset (Chr1[1..2 Mb], 5 Subclass_1 records) is deliberately below that
threshold, and their `short` dataset that does detect is 20 Mb of sequence - too
large to ship in a Tool Shed repository.

Measured thresholds, dante_tir 0.3.1 on MuDR/Mutator regions of tiny_pea:

    copies   fixture size   TIRs detected
       129        1.68 MB        0
       174        2.40 MB        4     <- what this script builds
       204        2.93 MB        5
       244        6.80 MB        7

A contiguous slice cannot reach 174 copies in a shippable size: the densest
MuDR/Mutator window in tiny_pea holds only 24 of them in 346 kb. So rather than
one slice, this takes the neighbourhood of each domain and concatenates them as
separate records - each becomes its own small "scaffold". That is ~8x more
copies per byte than a contiguous slice of the same size.

Regions are chosen densest-first, so REGIONS controls size and copy count
together. Every DANTE record that falls entirely inside a region is kept, not
just the MuDR/Mutator ones, so the input looks like a real DANTE annotation.
Coordinates are rewritten relative to each region.

The source dataset is in neither this repository nor upstream (81 MB FASTA,
22 MB GFF3). It lives at /mnt/ssd/dante_tir/test-data/ on the development
machine, so this script records how the fixtures were made rather than offering
to remake them from a clone.

Note that these tools must be run with their conda environment *activated*.
Invoking dante_tir.py by absolute path leaves cap3, mmseqs and blastn off PATH;
the run then still exits 0 and reports zero TIRs, and the real cause appears
only in <working_dir>/*.cap.err ("cap3: not found").

Usage:
    python3 make_test_fixture.py [SOURCE_DIR] [OUT_DIR]
"""
import collections
import sys
from pathlib import Path

SUPERFAMILY = "MuDR/Mutator"   # the most abundant Subclass_1 family in tiny_pea
REGIONS = 80                   # smallest tested count that still detects TIRs
FLANK = 7000                   # dante_tir extends ~6 kb, so a region needs more
MERGE_GAP = 20000              # copies closer than this share one region
WRAP = 60


def main():
    src = Path(sys.argv[1] if len(sys.argv) > 1 else "/mnt/ssd/dante_tir/test-data")
    out = Path(sys.argv[2] if len(sys.argv) > 2 else Path(__file__).parent)
    fasta_in, gff_in = src / "tiny_pea.fasta", src / "DANTE_tiny_pea.gff3"
    for p in (fasta_in, gff_in):
        if not p.exists():
            sys.exit(f"missing source file: {p}")

    # every DANTE record per chromosome, plus the seeds we build regions around
    by_chrom = collections.defaultdict(list)
    seeds = collections.defaultdict(list)
    with open(gff_in) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            s, e = int(f[3]), int(f[4])
            by_chrom[f[0]].append((s, e, f))
            if "Subclass_1" in f[8] and SUPERFAMILY in f[8]:
                seeds[f[0]].append((s, e))

    # merge nearby seeds, then keep the densest regions
    regions = []
    for c in seeds:
        cur = None
        for s, e in sorted(seeds[c]):
            if cur and s - cur[1] <= MERGE_GAP:
                cur = (cur[0], max(cur[1], e), cur[2] + 1)
            else:
                if cur:
                    regions.append((c, cur[0], cur[1], cur[2]))
                cur = (s, e, 1)
        if cur:
            regions.append((c, cur[0], cur[1], cur[2]))
    regions.sort(key=lambda r: (-r[3], r[0], r[1]))
    regions = sorted(regions[:REGIONS], key=lambda r: (r[0], r[1]))

    # load only the chromosomes the chosen regions need
    need = {r[0] for r in regions}
    seqs, name, buf = {}, None, []
    with open(fasta_in) as fh:
        for line in fh:
            if line.startswith(">"):
                if name in need:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name in need:
        seqs[name] = "".join(buf)

    kept = copies = bp = 0
    with open(out / "test_genome.fasta", "w") as fa, \
         open(out / "test_dante.gff3", "w") as gf:
        gf.write("##gff-version 3\n")
        for i, (c, s, e, _) in enumerate(regions, 1):
            lo, hi = max(1, s - FLANK), min(len(seqs[c]), e + FLANK)
            rid = f"{c}_r{i:03d}"
            sub = seqs[c][lo - 1:hi]
            fa.write(f">{rid}\n")
            for j in range(0, len(sub), WRAP):
                fa.write(sub[j:j + WRAP] + "\n")
            bp += len(sub)
            for rs, re_, f in by_chrom[c]:
                if rs < lo or re_ > hi:
                    continue
                g = list(f)
                g[0], g[3], g[4] = rid, str(rs - lo + 1), str(re_ - lo + 1)
                gf.write("\t".join(g) + "\n")
                kept += 1
                copies += "Subclass_1" in g[8] and SUPERFAMILY in g[8]

    print(f"{len(regions)} regions, {bp} bp ({bp / 1e6:.2f} Mb)")
    print(f"test_genome.fasta  {(out / 'test_genome.fasta').stat().st_size} bytes")
    print(f"test_dante.gff3    {(out / 'test_dante.gff3').stat().st_size} bytes, "
          f"{kept} DANTE records, {copies} {SUPERFAMILY} Subclass_1 copies")


if __name__ == "__main__":
    main()
