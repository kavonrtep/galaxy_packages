#!/usr/bin/env python3
"""Cut the DANTE test fixture out of the dante_ltr slice already in this repo.

The <tests> block of dante.xml referenced GEPY_test_long_1.fa and an expected
GFF3 that were never committed, so the test could not run from a clone - the
same fault as issue #5. Rather than add a new sequence, the fixture is a window
of dante_ltr/test-data/chr1_slice.fasta, so it is regenerable from this
repository alone.

The window 140,000-170,000 of that slice was chosen by counting DANTE records
per window in dante_ltr/test-data/dante_chr1_slice.gff3: it is the densest 30 kb
in the slice, with 18 complete protein domains forming two Ty3/gypsy
chromovirus/Tekay domain sets (GAG, PROT, RT, RH, INT, CHD). Six domain types
and one lineage are enough to assert on, and 30 kb keeps a LAST search against
REXdb to seconds.

Usage (from this directory):
    python3 make_test_fixture.py
"""
from pathlib import Path

HERE = Path(__file__).parent
SOURCE = HERE.parent.parent / "dante_ltr" / "test-data" / "chr1_slice.fasta"
START, END = 140_000, 170_000          # 0-based, half-open
WRAP = 60
OUT = HERE / "domains_slice.fasta"


def main():
    if not SOURCE.exists():
        raise SystemExit(f"missing source slice: {SOURCE}")
    seq = []
    with open(SOURCE) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    s = "".join(seq)[START:END]
    if len(s) != END - START:
        raise SystemExit(f"source is only {len(''.join(seq))} bp")
    with open(OUT, "w") as fh:
        fh.write(f">chr1_{START + 1}_{END} 30 kb of pea Chr1 with two Tekay domain sets\n")
        for i in range(0, len(s), WRAP):
            fh.write(s[i:i + WRAP] + "\n")
    print(f"{OUT.name}: {len(s)} bp, {OUT.stat().st_size} bytes")


if __name__ == "__main__":
    main()
