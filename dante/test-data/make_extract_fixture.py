#!/usr/bin/env python3
"""Build the input FASTA for the domains_extract tests (issue #6).

DANTE writes a '##No domains found' block in the middle of its GFF3 when one
part of the input has no hits, and dante_gff_to_dna.py of dante 0.2.10 skips
comment lines only at the top of the file. This fixture reproduces that with
three sequences, so the no-hit sequence sits between two with domains:

    seqA_domains  domains_slice.fasta (30 kb, 18 domains)
    seqB_nohits   20 kb of fixed-seed random sequence (no hits)
    seqC_domains  dante_ltr/test-data/chr1_slice.fasta 60,000-90,000 (3 domains)

The GFF3 fixtures are DANTE's own output on that FASTA, from dante 0.2.10:

    dante --query extract_input.fasta --domain_gff extract_dante_full.gff3 \\
          --database Viridiplantae_v4.0 --scoring_matrix BL80 --cpu 3
    grep '^#' extract_dante_full.gff3 | head -3 > extract_header_only.gff3

With three CPUs DANTE runs the sequences as separate parts, which is what puts
the comment block mid-file (lines 25-28).

Usage (from this directory):
    python3 make_extract_fixture.py
"""
import random
from pathlib import Path

HERE = Path(__file__).parent
WRAP = 60


def read_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))


def main():
    a = read_seq(HERE / "domains_slice.fasta")
    c = read_seq(HERE.parent.parent / "dante_ltr" / "test-data" / "chr1_slice.fasta")[60_000:90_000]
    rng = random.Random(6)
    b = "".join(rng.choice("ACGT") for _ in range(20_000))
    with open(HERE / "extract_input.fasta", "w") as out:
        for name, seq in (("seqA_domains", a), ("seqB_nohits", b), ("seqC_domains", c)):
            out.write(f">{name}\n")
            for i in range(0, len(seq), WRAP):
                out.write(seq[i:i + WRAP] + "\n")


if __name__ == "__main__":
    main()
