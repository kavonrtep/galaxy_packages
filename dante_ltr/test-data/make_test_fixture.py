#!/usr/bin/env python3
"""Build the dante_ltr test fixtures by slicing tiny_pea.

A dante_ltr fixture needs two things that pull in opposite directions: enough
sequence for complete LTR retrotransposons to be recognised by their structure,
and small enough to ship in a Tool Shed repository. Measured with dante_ltr
0.6.3.1 on contiguous slices of tiny_pea Chr1, lineage mode:

    slice     DANTE records   elements detected   fixture size
      500 kb             58                  12        0.56 MB   <- this script
     1000 kb             96                  18        1.10 MB
     1500 kb            201                  39        1.70 MB
     2000 kb            400                  80        2.37 MB

500 kb already yields 12 elements, so there is no reason to ship more. All four
tools work on it: clean_ltr.R produces all six outputs, dante_ltr_summary renders
its HTML, and dante_ltr_to_library returns 11 representative sequences even at
--min_coverage 3. Core mode finds the same 12 elements as lineage mode.

**The slice must be contiguous.** The dante_tir fixture in this repository is 80
concatenated regions, which is far more copies per byte, but clean_ltr.R dies on
that shape in multiplicity_of_te with "non-numeric argument to binary operator" -
it computes coverage across the whole sequence. See dante_tir/test-data/README.md.

Three files are produced. The third is the output of the first two, kept so the
three downstream tools have a realistic input without having to chain tools in a
test:

    chr1_slice.fasta            the genome slice, one record named Chr1
    dante_chr1_slice.gff3       every DANTE record inside it
    dante_ltr_chr1_slice.gff3   dante_ltr's own output for those two

The source dataset is in neither this repository nor upstream (81 MB FASTA,
22 MB GFF3). It lives at /mnt/ssd/dante_tir/test-data/ on the development
machine, so this script records how the fixtures were made rather than offering
to remake them from a clone.

Usage:
    python3 make_test_fixture.py [SOURCE_DIR] [OUT_DIR]

Then regenerate the third file by running dante_ltr on the first two:
    dante_ltr -g dante_chr1_slice.gff3 -s chr1_slice.fasta -o out -c 4
    cp out.gff3 dante_ltr_chr1_slice.gff3
"""
import sys
from pathlib import Path

CHROM = "Chr1"
END = 500_000          # smallest tested slice that still detects elements
WRAP = 60


def main():
    src = Path(sys.argv[1] if len(sys.argv) > 1 else "/mnt/ssd/dante_tir/test-data")
    out = Path(sys.argv[2] if len(sys.argv) > 2 else Path(__file__).parent)
    fasta_in, gff_in = src / "tiny_pea.fasta", src / "DANTE_tiny_pea.gff3"
    for p in (fasta_in, gff_in):
        if not p.exists():
            sys.exit(f"missing source file: {p}")

    seq, wanted = [], False
    with open(fasta_in) as fh:
        for line in fh:
            if line.startswith(">"):
                if wanted:
                    break
                wanted = line[1:].split()[0] == CHROM
                continue
            if wanted:
                seq.append(line.strip())
    s = "".join(seq)[:END]
    if len(s) < END:
        sys.exit(f"{CHROM} is shorter than {END} bp")

    with open(out / "chr1_slice.fasta", "w") as fh:
        fh.write(f">{CHROM}\n")
        for i in range(0, len(s), WRAP):
            fh.write(s[i:i + WRAP] + "\n")

    kept = ltr = 0
    with open(gff_in) as fin, open(out / "dante_chr1_slice.gff3", "w") as fout:
        fout.write("##gff-version 3\n")
        for line in fin:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != CHROM or int(f[4]) > END:
                continue
            fout.write(line)
            kept += 1
            ltr += "Class_I|LTR" in f[8]

    print(f"{CHROM}:1-{END}  {len(s)} bp")
    print(f"chr1_slice.fasta       {(out / 'chr1_slice.fasta').stat().st_size} bytes")
    print(f"dante_chr1_slice.gff3  {kept} DANTE records, {ltr} of them Class_I|LTR")


if __name__ == "__main__":
    main()
