#!/usr/bin/env python3
"""Generate a deterministic synthetic genome with one clean high-copy tandem
array, so TideCluster reliably produces TRC_1 and a non-empty consensus dimer
library. Fixed seed => identical bytes on every regeneration."""
import random

random.seed(1234)  # deterministic
BASES = "ACGT"

def rand_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))

def mutate(seq, rate):
    out = []
    for b in seq:
        if random.random() < rate:
            out.append(random.choice(BASES))   # substitution
        else:
            out.append(b)
    return "".join(out)

# One satellite: 172 bp monomer, ~500 copies, ~3% per-base divergence per copy (array ~86 kb, above the default
monomer = rand_seq(172)
array = "".join(mutate(monomer, 0.03) for _ in range(500))   # min_total_length=50000 TAREAN threshold)

# Unique random flanks so the array sits inside a realistic contig.
left = rand_seq(25000)
right = rand_seq(25000)
seq = left + array + right

def wrap(s, w=60):
    return "\n".join(s[i:i+w] for i in range(0, len(s), w))

with open("synthetic_satellite.fasta", "w") as f:
    f.write(">chr_test synthetic contig with one 172bp satellite array\n")
    f.write(wrap(seq) + "\n")

print(f"monomer={len(monomer)}bp copies=500 array={len(array)}bp total={len(seq)}bp")
