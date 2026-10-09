#!/usr/bin/env python3
"""Build the tidecluster_annotation fixtures out of comparative_run_a.zip.

The annotation step needs something the comparative fixtures do not provide on
their own: a run directory carrying the prefix the wrapper passes (-pr
tidecluster) and a `<prefix>_consensus/TRC*dimers.fasta` set, which is the input
TideCluster.py annotation globs for. Its own clustering GFF3 has no
`consensus_sequence` attributes, so the fallback path in TideCluster.py is not
available either.

Rather than ship a second copy of a real run, this derives both inputs from
comparative_run_a.zip, which is already in this directory:

    annotation_run.zip          tidecluster_clustering.gff3 (6 regions) and
                                tidecluster_consensus/TRC_{26,28,36}_dimers.fasta
                                (9 sequences, 11 kb of sequence)
    annotation_library.fasta    two references, each the first dimer of one of
                                those clusters, relabelled in RepeatMasker format

Three clusters, two of them represented in the library, is the smallest shape
that covers the three outcomes the tool has: a cluster annotated from the
library, a second one annotated as a different class, and a third with no
reference at all, which has to come back as `annotation=NA` rather than as a
guess. The classes are deliberately made up (`Satellite/testfam`, `rDNA/45S`)
and taken only from the FASTA headers, which is what pins the documented
behaviour that the string after `#` is what gets reported.

The library references are per-array dimers rather than TAREAN consensus
sequences, because the dimer library in comparative_run_a.zip only covers
TRC_1, 4, 6, 8 and 11 - all of them larger clusters. A reference taken from a
cluster's own pool is a near-exact match for its siblings, which is what makes
the expected annotation deterministic.

Usage (from this directory):
    python3 make_annotation_fixture.py
"""
import re
import shutil
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).parent
SOURCE = HERE / "comparative_run_a.zip"
KEEP = ["26", "28", "36"]
# cluster -> header for the library reference built from its first dimer
REFS = [("26", "testsat_A#Satellite/testfam"), ("28", "testsat_B#rDNA/45S")]
# A fixed timestamp so regenerating gives byte-identical archives.
ZIP_DATE = (1980, 1, 1, 0, 0, 0)
WRAP = 60


def read_fasta(path):
    recs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), WRAP):
                fh.write(seq[i:i + WRAP] + "\n")


def main():
    if not SOURCE.exists():
        sys.exit(f"missing source archive: {SOURCE}")
    work = HERE / "_annotation_fixture_tmp"
    shutil.rmtree(work, ignore_errors=True)
    src, run = work / "src", work / "run"
    (run / "tidecluster_consensus").mkdir(parents=True)
    src.mkdir(parents=True)
    with zipfile.ZipFile(SOURCE) as z:
        z.extractall(src)

    pool = read_fasta(src / "tc_consensus" / "consensus_sequences_all.fasta")
    for trc in KEEP:
        sel = [(n, s) for n, s in pool if n.split("_")[1] == trc]
        if not sel:
            sys.exit(f"TRC_{trc} not present in the source consensus pool")
        write_fasta(run / "tidecluster_consensus" / f"TRC_{trc}_dimers.fasta", sel)
        print(f"TRC_{trc}: {len(sel)} sequences, {sum(len(s) for _, s in sel)} bp")

    kept = 0
    with open(run / "tidecluster_clustering.gff3", "w") as out:
        out.write("##gff-version 3\n")
        with open(src / "tc_clustering.gff3") as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                m = re.search(r"Name=TRC_(\d+)", line)
                if m and m.group(1) in KEEP:
                    out.write(line)
                    kept += 1
    print(f"tidecluster_clustering.gff3: {kept} regions")

    library = []
    for trc, header in REFS:
        name, seq = read_fasta(
            run / "tidecluster_consensus" / f"TRC_{trc}_dimers.fasta")[0]
        library.append((header, seq))
        print(f"library: {header} <- {name} ({len(seq)} bp)")
    write_fasta(HERE / "annotation_library.fasta", library)

    archive = HERE / "annotation_run.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
        for path in sorted(p for p in run.rglob("*") if p.is_file()):
            info = zipfile.ZipInfo(str(path.relative_to(run)), date_time=ZIP_DATE)
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, path.read_bytes())
    shutil.rmtree(work)
    print(f"{archive.name}: {archive.stat().st_size} bytes")


if __name__ == "__main__":
    main()
