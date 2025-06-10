#!/usr/bin/env python3
import argparse
import sys
from argparse import ArgumentDefaultsHelpFormatter
from collections import namedtuple
from collections import OrderedDict
from Bio import SeqIO

SEQ_FORMAT = "fasta"
# Default sampling, used on argparse:
DEFAULT_READ_LENGTH = 200
DEFAULT_INSERT_LENGTH = 700
DEFAULT_COVERAGE = 0.1

Sequences_summary = namedtuple('Fasta_summary',
                               ['total_length', 'number_of_sequence',
                                'id_length', 'file_path', 'format'])

Coordinates = namedtuple('Coordinates', "id start1 end1 start2 end2")


def get_sequences_summary(seq_file):
    ''' return basic characteristic of sequences '''
    id_length = OrderedDict()
    totat_length = 0
    N = 0
    for seqs in SeqIO.parse(seq_file, SEQ_FORMAT):
        id_length[seqs.id] = len(seqs)
        totat_length += len(seqs)
        N += 1
    return Sequences_summary(totat_length, N, id_length, seq_file, SEQ_FORMAT)


def get_short_pseudoreads_position(fasta_summary, sampling_options):
    """Return selected position on long read sequences
    Arguments:
    fasta_summary - namedtuple Fasta_summaty containing information about sequences
    sampling options - namedtuple, specified how sequences should be sampled
    Return value:
    (sequence_id, start1, end1, start2, end2)
    """
    interval = int(2 * sampling_options.read_length /
                   sampling_options.coverage)
    for seqname, length in fasta_summary.id_length.items():
        start_positions = range(1, length, interval)
        for s in start_positions:
            yield Coordinates(seqname, s, s + sampling_options.read_length,
                              s + sampling_options.insert_length -
                              sampling_options.read_length,
                              s + sampling_options.insert_length)


def extract_short_reads(summary, args):
    '''yield short reades sampled from long reads
    Arguments:
    summary.. named tuple specifie sequences properties, path, length, idslist
    args ..... Define how short sequences should be generated
    '''
    pos = get_short_pseudoreads_position(summary, args)
    coords = next(pos)
    index = 0
    for i in SeqIO.parse(summary.file_path, summary.format):
        index += 1
        while True:
            if coords.id == i.id:
                # forward read
                subseq_f = i[coords.start1:coords.end1]
                subseq_f.id = "{}_{}_{}_f".format(index, coords.start1,
                                                  coords.end1)
                subseq_f.description = ""
                # reverse complement read
                subseq_r = i[coords.start2:coords.end2].reverse_complement()
                subseq_r.id = "{}_{}_{}_r".format(index, coords.start1,
                                                  coords.end1)
                subseq_r.description = ""
                # return only if sequences are long enough
                if len(subseq_r) == args.read_length:
                    yield subseq_f
                    yield subseq_r
                coords = next(pos)
            else:
                break


def long2short(args):
    '''Sample short reads from long sequences
    args contain these attributes::
    ------------
    input_file   - path to file in fasta format
    output_file  - path to output file, fasta format
    options      - options is named tuple and specifies read length
                   coverage, insert length, max number of sequences which will be return

    '''
    summary = get_sequences_summary(args.input.name)
    with open(args.output.name, 'w') as f:
        for i in extract_short_reads(summary, args):
            SeqIO.write(i, f, SEQ_FORMAT)


def get_args():
    '''Parses command line arguments '''
    description = "Creates pseudo short reads from long oxford nanopore reads"
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i',
                        '--input',
                        type=argparse.FileType('r'),
                        help="file with long reads in fasta format")
    parser.add_argument('-o',
                        '--output',
                        type=argparse.FileType('w'),
                        help="Output file name")
    parser.add_argument("-cov",
                        "--coverage",
                        type=float,
                        default=DEFAULT_COVERAGE,
                        help="samplig coverage")
    parser.add_argument(
        "-L",
        "--insert_length",
        type=int,
        default=DEFAULT_INSERT_LENGTH,
        help="length of insert, must be longer than read length")
    parser.add_argument("-l",
                        "--read_length",
                        type=int,
                        default=DEFAULT_READ_LENGTH,
                        help="read length")

    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    #hassert args.insert_length > args.read_length, "read length must be shorter than insert length"
    return args


def main():
    '''Sample short reads from long sequences
    Files path are passed as command line positional arguments
                        '''
    args = get_args()
    long2short(args)


if __name__ == "__main__":
    main()
