#!/usr/bin/env python3
"""Cut the FASTQ fixtures for the two read-filtering tools.

The <tests> blocks of paired_fastq_filtering.xml and single_fastq_filtering.xml
named ERR215189_1_part.fastq.gz and friends, none of which were ever committed,
so neither test could run from a clone - the same fault as issue #5.

The pair here is the first 1000 read pairs of SRR089356_{1,2}.fastq.gz, which is
already in this directory but is excluded from the published tarball at 7.8 MB
per side. 1000 pairs of 100 bp reads come to about 80 kB per side, small enough
to ship, and enough for the tools to do real work: at the default quality
settings part of the input is removed, so the tests see filtering rather than a
pass-through.

Taking the first N records rather than a random sample keeps the fixture
reproducible without carrying a seed, and gzip is written with mtime 0 so
regenerating gives byte-identical files.

Usage (from this directory):
    python3 make_filtering_fixture.py
"""
import gzip
from pathlib import Path

HERE = Path(__file__).parent
PAIRS = 1000
SOURCES = {"1": HERE / "SRR089356_1.fastq.gz", "2": HERE / "SRR089356_2.fastq.gz"}
OUT = "filtering_{side}.fastq.gz"


def main():
    for side, src in SOURCES.items():
        if not src.exists():
            raise SystemExit(f"missing source: {src}")
        lines = []
        with gzip.open(src, "rt") as fh:
            for i, line in enumerate(fh):
                if i >= PAIRS * 4:
                    break
                lines.append(line)
        if len(lines) < PAIRS * 4:
            raise SystemExit(f"{src.name} holds fewer than {PAIRS} reads")
        out = HERE / OUT.format(side=side)
        # mtime=0 so the bytes do not change between regenerations
        with gzip.GzipFile(out, "wb", mtime=0) as fh:
            fh.write("".join(lines).encode())
        print(f"{out.name}: {PAIRS} reads, {out.stat().st_size} bytes")


if __name__ == "__main__":
    main()
