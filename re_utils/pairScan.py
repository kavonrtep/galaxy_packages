#!/usr/bin/env python
import sys
import os
from optparse import OptionParser
import Levenshtein


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


def readSingleSeq(file):
    line = file.readline()
    if not line:
        return False  # end of file
    if line[0] != ">":
        raise Error("no header on the first line")
    seqname = line[1:].strip()
    seq = ""
    # read sequences
    while True:
        last_pos = file.tell()
        line = file.readline()
        if not line:
            break
        if line[0] == ">":
            file.seek(last_pos)
            break
        seq = seq + line.strip()
    return {'name': seqname, 'sequence': seq}


def writeSingleSeq(fileobject, seq):
    fileobject.write(">")
    fileobject.write(seq['name'] + "\n")
    fileobject.write(seq['sequence'] + "\n")


def comparePairs(seq1, seq2, max_mismatch=3, offset=5):
    s1 = seq1['sequence'].lower()
    s2 = seq2['sequence'].lower()[::-1]
    m = 0
    intab = "ctagn"
    outtab = "gatcn"
    trantab = str.maketrans(intab, outtab)
    s2 = s2.translate(trantab)
    s1 = "-" * offset + s1
    s2 = s2 + "-" * offset
    n1 = len(s1)
    n2 = len(s2)
    m = 0
    for i in range(1, min(n1 + 1, n2 + 1)):
        #remove tails is any:
        ss1 = s1[n1 - i:n1]
        ss2 = s2[0:i]
        added = ss1.count("-") + ss2.count("-")
        d = Levenshtein.hamming(ss1, ss2) - added
        if 100.0 * d / i <= max_mismatch:
            m = max(m, i - d - added)
    return m


def split_file(filename, N, min_chunk=2):
    f1 = open(filename, 'r')
    filenames = [filename + "." + str(i) for i in range(N)]
    f2 = list(map(open, filenames, 'w' * N))
    while True:
        for i in f2:
            for j in range(min_chunk):
                line = f1.readline()
                if not line:
                    [i.close() for i in f2]
                    f1.close()
                    return filenames
                i.write(line)


def find_overlapping_sequences(seqfile,
                               seqfile2=None,
                               seqfile_good="",
                               seqfile_bad="",
                               min_overlap=30,
                               max_mismatch=2,
                               offset=5):
    ''' return id ove overlaping pairs - only first id is returned '''
    # default names - if empty
    if seqfile_good == "":
        seqfile_good = seqfile + ".pass"
    if seqfile_bad == "":
        seqfile_bad = seqfile + ".bad"

    minscore = min_overlap * 2

    fgood = open(seqfile_good, 'w')
    fbad = open(seqfile_bad, 'w')
    f = open(seqfile, 'r')
    if seqfile2:
        f2 = open(seqfile2)
    else:
        f2 = f
    while True:
        seq1 = readSingleSeq(f)
        seq2 = readSingleSeq(f2)
        if not seq1 or not seq2:
            break  # end of file
        score = comparePairs(seq1, seq2, max_mismatch, offset=offset)
        if score > min_overlap:
            writeSingleSeq(fbad, seq1)
            writeSingleSeq(fbad, seq2)
        else:
            writeSingleSeq(fgood, seq1)
            writeSingleSeq(fgood, seq2)
    f.close()
    if not f2.closed:
        f2.close
    fgood.close()
    fbad.close()


def main():
    parser = OptionParser()
    parser.add_option("-f",
                      "--fasta_file",
                      dest="seqfile",
                      help="input sequences in fasta format")
    parser.add_option(
        "-r",
        "--fasta_file2",
        default=None,
        dest="seqfile2",
        help=
        "input sequences in fasta format, second file should be specified if pairs are not interlaced, all pairs must be complete!")
    parser.add_option("-p",
                      "--fasta_file_pass",
                      dest="seqfile_good",
                      help="output file with good sequences",
                      default='')
    parser.add_option("-b",
                      "--fasta_file_bad",
                      dest="seqfile_bad",
                      help="output file with bad sequences",
                      default='')
    parser.add_option("-o",
                      "--minimal_overlap",
                      dest="min_overlap",
                      help="minimal overlap between pair ends",
                      default='30')
    parser.add_option(
        "-m",
        "--max_mismatch",
        dest="max_mismatch",
        help="maximum number of mismatches in overlap per 100 nt",
        default='2')
    parser.add_option("-s",
                      "--offset",
                      dest="offset",
                      help="maximum offset",
                      default='5')
    options, args = parser.parse_args()
    find_overlapping_sequences(seqfile=options.seqfile,
                               seqfile2=options.seqfile2,
                               seqfile_good=options.seqfile_good,
                               seqfile_bad=options.seqfile_bad,
                               min_overlap=int(options.min_overlap),
                               max_mismatch=int(options.max_mismatch),
                               offset=int(options.offset))


if __name__ == "__main__":
    main()
