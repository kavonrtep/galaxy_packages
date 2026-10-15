#!/usr/bin/env python3
"""Build the two fixtures for the carp input-check tests from genome_micro.fasta.

input_check_reads.fasta       600 reads of 150 bp sampled from the genome
input_check_fragmented.fasta  the genome cut into 1.5-4.5 kb pieces (N50 < 5 kb)

Fixed seed, so the output is the same on every run.
"""
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path, records, width=60):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def main():
    rng = random.Random(20261013)
    genome = read_fasta(os.path.join(HERE, "genome_micro.fasta"))
    names = sorted(genome)

    reads = []
    for i in range(600):
        seq = genome[rng.choice(names)]
        start = rng.randrange(len(seq) - 150)
        reads.append((f"read_{i + 1}", seq[start:start + 150]))
    write_fasta(os.path.join(HERE, "input_check_reads.fasta"), reads, width=150)

    pieces = []
    for name in names:
        seq, pos = genome[name], 0
        while pos < len(seq):
            size = rng.randint(1500, 4500)
            pieces.append((f"{name}_ctg{len(pieces) + 1}", seq[pos:pos + size]))
            pos += size
    write_fasta(os.path.join(HERE, "input_check_fragmented.fasta"), pieces)


if __name__ == "__main__":
    main()
