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

# Reference library for the tc_reannotate ("Annotate Genome") test: the
# unmutated monomer, named in RepeatMasker format so the class after "#" is
# what the tool reports. Written from the monomer computed above rather than
# from a copy taken out of the array, and after the genome so the random stream
# - and therefore synthetic_satellite.fasta - is unchanged by its presence.
with open("trc_library.fasta", "w") as f:
    f.write(">TRC_1#Satellite/synthetic\n")
    f.write(wrap(monomer) + "\n")

print(f"monomer={len(monomer)}bp copies=500 array={len(array)}bp total={len(seq)}bp")
print("trc_library.fasta: 1 reference, TRC_1#Satellite/synthetic")
