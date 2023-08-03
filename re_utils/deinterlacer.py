#!/usr/bin/env python3
'''very simple deinterlacer - fasta and fastq'''
import sys
import itertools


def is_header(line, counter, fasta):
    ''' return True is line is header '''
    if fasta:
        if line[0] == ">":
            return True
    else:
        if counter == 4 and line[0] == "@":
            return True
    return False


def main():
    '''deinterlace fasta or fastq format'''
    infile = sys.argv[1]
    file_a = sys.argv[2]
    file_b = sys.argv[3]
    with open(infile) as f, open(file_a, 'w') as A, open(file_b, 'w') as B:
        ABiter = itertools.cycle([A, B])
        counter = 3  # four lines per record in fastq
        pos = f.tell()
        is_fasta = f.readline()[0] == ">"
        f.seek(pos)
        for line in f:
            counter += 1
            if is_header(line, counter, is_fasta):
                fout = next(ABiter)
                counter = 0
            fout.write(line)


if __name__ == "__main__":
    main()
