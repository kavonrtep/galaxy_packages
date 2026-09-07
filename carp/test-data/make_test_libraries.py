#!/usr/bin/env python3
"""Generate the two optional repeat libraries used by the CARP tool test.

Both files exist to exercise the classification validation CARP 1.7.1 added:
`custom_library` and `tandem_repeat_library` are the only classification-bearing
inputs the pipeline does not produce itself, and a header whose class is not in
CARP's vocabulary aborts the run at `validate_classifications`. The fixtures
therefore carry canonical classes only, so the test fails if that contract is
broken in either direction.

Deterministic: fixed seed, no dependencies. Regenerate with

    python3 make_test_libraries.py

Verify the headers against the pipeline's own vocabulary with

    classification.py validate --mode fasta --source RepeatMasker <file>
"""

import random

SEED = 20260907
LINE = 60

# Canonical classes, as accepted by classification_vocabulary.yaml. The three
# TE classes come from `classification.py list-canonical`; `Satellite` is a
# RepeatMasker-native class the validator also accepts and is the natural label
# for a tandem-repeat library.
CUSTOM = [
    ("myAle_1", "Class_I/LTR/Ty1_copia/Ale", 600),
    ("myCACTA_1", "Class_II/Subclass_1/TIR/EnSpm_CACTA", 500),
    ("my18S_1", "rDNA/45S_rDNA/18S", 400),
]
TANDEM = [
    ("TestTR_A", "Satellite", 55, 8),
    ("TestTR_B", "Satellite", 90, 5),
]


def wrap(seq):
    return "\n".join(seq[i:i + LINE] for i in range(0, len(seq), LINE))


def random_dna(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def main():
    rng = random.Random(SEED)

    with open("custom_library.fasta", "w") as fh:
        for name, cls, length in CUSTOM:
            fh.write(f">{name}#{cls}\n{wrap(random_dna(rng, length))}\n")

    # Satellites: a monomer tiled head-to-tail, which is what a tandem-repeat
    # reference library actually holds.
    with open("tandem_repeat_library.fasta", "w") as fh:
        for name, cls, monomer_len, copies in TANDEM:
            monomer = random_dna(rng, monomer_len)
            fh.write(f">{name}#{cls}\n{wrap(monomer * copies)}\n")


if __name__ == "__main__":
    main()
