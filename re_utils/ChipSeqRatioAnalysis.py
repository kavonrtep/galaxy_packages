#!/usr/bin/env python3
import re
import argparse
import csv
import os
import sys
import os.path
import subprocess
import shlex
import multiprocessing as mp
import tempfile
import itertools as it


def get_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument('-m',
                        '--max_cl',
                        default=300,
                        type=int,
                        help='Sets the maximum cluster number. Default = 300')
    parser.add_argument('-b',
                        '--bitscore',
                        default=30,
                        type=int,
                        help='minimal bitscore to report')
    parser.add_argument(
        '-n',
        '--nproc',
        default=mp.cpu_count(),
        type=int,
        help='Sets the number of cpus to be used. Default = all available')
    parser.add_argument('-c',
                        '--ChipSeq',
                        required=True,
                        help='Fasta file containing the Chip Sequences')
    parser.add_argument('-i',
                        '--InputSeq',
                        required=True,
                        help='Fasta file containing the Input Sequences')
    parser.add_argument(
        '-o',
        '--output',
        required=False,
        default='ChipSeqRatio.csv',
        help=('Specify a name for the CSV file to which the'
              ' output will be save to. Default: ChipSeqRatio.csv'))
    parser.add_argument(
        '-ht',
        '--html',
        required=False,
        default='ChipSeqRatioReport',
        help='Specify a name for the html report. Default : ChipSeqRatioReport')
    parser.add_argument('-k',
                        '--Contigs',
                        required=True,
                        help='Contig file for blast')
    args = parser.parse_args()
    return args


def split(filename, num_chunks):
    # splits a files into nproc files
    files = []
    temp_files = []

    # creating a list with nproc temporary files
    for f in range(num_chunks):
        temp_files.append(tempfile.NamedTemporaryFile('w', delete=False))
    with open(filename, 'r') as f:
        end = os.path.getsize(filename)
        for temp in it.cycle(temp_files):
            # cycling indefenitly through temp files
            ID = f.readline()  # get sequence id
            if ID[0] == '>':
                temp.write(ID)
                SEQP = f.tell()  # get file pointer location
                while f.readline()[0] is not '>':
                    f.seek(SEQP)  # jump to last saved file pointer
                    temp.write(f.readline())  # write sequen
                    SEQP = f.tell()  # overwrite last file pointer location
                    # break loop if file pointer reached the EOF
                    if (SEQP == end):
                        break
            if (SEQP == end):  # break loop if file pointer reached the EOF
                break
            f.seek(SEQP)

    for f in range(num_chunks):
        # save temp files names into list for further use
        files.append(temp_files[f].name)
        temp_files[f].close()  # close temp files
    return files


def blast(query, database, bitscore):
    # blast a file for given arguments and save result to a list
    print(query)
    arguments = (
        "blastn -task blastn -db {} -query {} "
        "-evalue 1e-2 -gapopen 5 -gapextend 2 -word_size 11 -num_alignments 1"
        " -penalty -3 -reward 2 -outfmt 6 -dust no").format(database, query)
    cmd = shlex.split(arguments)
    Blast_Output = [0 for x in range(max_cl + 1)]
    ma = re.compile('(\S+)\tCL(\d+)Contig')  # expression to check for
    with subprocess.Popen(cmd,
                          stdout=subprocess.PIPE,
                          universal_newlines=True) as p:
        for line in p.stdout:
            if float(line.split()[11]) > bitscore:
                gr = ma.match(line)
                previous_query = ''
                if gr:
                    if (gr.group(1) != previous_query):
                        if (int(gr.group(2)) > max_cl):
                            Blast_Output[0] = Blast_Output[0] + 1
                        else:
                            Blast_Output[int(gr.group(2))] = Blast_Output[
                                int(gr.group(2))] + 1
                        previous_query = gr.group(1)
    return Blast_Output


def ReduceLists(x, y):
    ''' reduces two lists into a 2-dim matrix '''
    Matrix = [[0 for i in range(max_cl + 1)] for i in range(2)]
    for i in range(len(x)):
        for j in range(len(x[i])):
            Matrix[0][j] = Matrix[0][j] + x[i][j]
    for i in range(len(y)):
        for j in range(len(y[i])):
            Matrix[1][j] = Matrix[1][j] + y[i][j]
    return Matrix


def fasta_size(fastafile):
    with open(fastafile, 'r') as f:
        s = 0
        for i in f:
            if i[0] == ">":
                s += 1
    return s


def makeblastdb(filename):
    dbtmp = tempfile.NamedTemporaryFile()
    cmd = [
        'makeblastdb', '-in', filename, '-input_type', 'fasta', '-dbtype',
        'nucl', '-out', dbtmp.name
    ]
    subprocess.call(cmd)
    return dbtmp


if __name__ == "__main__":
    args = get_arguments()
    max_cl = args.max_cl
    output = args.output
    HTMLreport = args.html
    contigs = args.Contigs

    # Creation of database
    db = makeblastdb(contigs)

    inputN = fasta_size(args.InputSeq)
    chipN = fasta_size(args.ChipSeq)

    # Reading and distribution of data to temp files for multiprocessing
    filesC = split(args.ChipSeq, args.nproc)
    filesI = split(args.InputSeq, args.nproc)

    # start of parallized blast
    pool = mp.Pool(processes=args.nproc)
    results = [pool.apply_async(blast, args=(f, db.name, args.bitscore)) for f in filesC]
    Cout = [p.get() for p in results]
    results = [pool.apply_async(blast, args=(f, db.name, args.bitscore)) for f in filesI]
    Iout = [p.get() for p in results]

    # Merging of blast output into a 2-dim matrix
    Matrix = ReduceLists(Cout, Iout)

    with open(args.output, 'w') as f:
        print("Cluster", "Chip_Hits", "Input_Hits", sep='\t', file=f)
        for hit in range(1, args.max_cl + 1):
            print(hit, Matrix[0][hit], Matrix[1][hit], sep='\t', file=f)
    Rarguments = "Rscript " + \
        os.path.dirname(__file__) + "/ChipSeqRatioAnalysis.R"
    # order is important - programmed by georg - this it realy ugly!
    args = shlex.split(Rarguments)
    args.append(output)
    args.append(HTMLreport)
    args.append(str(inputN))
    args.append(str(chipN))
    with subprocess.Popen(args, stderr=subprocess.PIPE) as p:
        print("Creating HTML report")
        stdout, stderr = p.communicate()
        if (len(stderr) > 0):
            print(stderr)
    # cleanup
    for i in filesC + filesI:
        os.unlink(i)
