#!/usr/bin/env python3
"""Check that a FASTA file is a genome assembly CARP can annotate.

Run before the pipeline, so that unsuitable input fails in about a minute
instead of hours into a run. Two cases are distinguished:

- sequencing reads (wrong input): always an error;
- a very fragmented assembly (N50 below --min-n50): an error unless
  --allow-fragmented is given, in which case a warning is printed.

Statistics go to stdout, messages to stderr. Exit status 1 on rejection.
"""
import argparse
import collections
import sys

# Reads: no sequence longer than this, or at least READ_MODE_FRAC of the
# sequences share one length no longer than this. Short reads are at most
# 300 bp, about 600 bp for merged pairs.
READ_MAX_LEN = 1000
READ_MODE_FRAC = 0.5


def sequence_lengths(path):
    lengths = []
    current = None
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if current is not None:
                    lengths.append(current)
                current = 0
            elif current is not None:
                current += len(line.strip())
    if current is not None:
        lengths.append(current)
    return lengths


def n50(lengths, total):
    acc = 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if 2 * acc >= total:
            return length
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fasta")
    parser.add_argument("--min-n50", type=int, default=5000,
                        help="minimum N50 in bp (default 5000)")
    parser.add_argument("--allow-fragmented", action="store_true",
                        help="warn instead of failing when N50 is below --min-n50")
    args = parser.parse_args()

    lengths = sequence_lengths(args.fasta)
    if not lengths:
        sys.exit("CARP input check: the genome file contains no FASTA sequences.")

    n = len(lengths)
    total = sum(lengths)
    longest = max(lengths)
    mode_len, mode_count = collections.Counter(lengths).most_common(1)[0]
    mode_frac = mode_count / n
    asm_n50 = n50(lengths, total)

    stats = (
        f"  sequences:           {n:,}\n"
        f"  total length:        {total:,} bp\n"
        f"  longest sequence:    {longest:,} bp\n"
        f"  N50:                 {asm_n50:,} bp\n"
        f"  most common length:  {mode_len:,} bp ({mode_frac:.1%} of sequences)\n"
    )
    if longest <= READ_MAX_LEN or (mode_len <= READ_MAX_LEN and mode_frac >= READ_MODE_FRAC):
        sys.stderr.write(
            "CARP input check FAILED: the input looks like sequencing reads, not a "
            "genome assembly.\n" + stats +
            f"No sequence is longer than {READ_MAX_LEN:,} bp, or most sequences share "
            f"one short length. CARP annotates genome assemblies; assemble the reads "
            f"first and run CARP on the assembly.\n")
        return 1

    if asm_n50 < args.min_n50:
        detail = (
            f"The assembly N50 is {asm_n50:,} bp, below {args.min_n50:,} bp. CARP "
            f"identifies full-length transposable elements and long tandem repeat "
            f"arrays, which do not fit in contigs this short: expect few complete "
            f"elements, a poor repeat library, and possibly a failed run.\n")
        if args.allow_fragmented:
            sys.stderr.write(
                "CARP input check WARNING: very fragmented assembly; continuing "
                "because the override is set.\n" + stats + detail)
            return 0
        sys.stderr.write(
            "CARP input check FAILED: the assembly is too fragmented for CARP.\n"
            + stats + detail +
            "To run anyway, enable 'Run on a very fragmented assembly' under "
            "Advanced settings.\n")
        return 1

    print("CARP input check passed:\n" + stats, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
