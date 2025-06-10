#!/usr/bin/env python3
import sys
import argparse
import random
from argparse import ArgumentDefaultsHelpFormatter
from Bio import SeqIO
from long2short import get_sequences_summary


def get_args():
    '''Parses command line arguments '''
    description = ("Create sample of long reads, instead"
                   " of setting number of reads to be sampled,"
                   "total length of all sampled sequences is defined")

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
    parser.add_argument("-l",
                        "--total_length",
                        type=int,
                        help="total length of sampled output")

    parser.add_argument("-s",
                        "--seed",
                        type=int,
                        default=123,
                        help="random number generator seed")

    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return args


def make_sequence_sample(args, seq_summary):
    '''
    create output sequence of required length or larger, save to fasta
    '''
    random.seed(args.seed)
    print(args.seed)
    selected_seq_id = set()
    cummulative_length = 0
    for seq_id in random.sample(seq_summary.id_length.keys(),
                                seq_summary.number_of_sequence):
        selected_seq_id.add(seq_id)
        cummulative_length += seq_summary.id_length[seq_id]
        if cummulative_length >= args.total_length:
            break
    # to make it efficient same orger is necesary
    selected_seq_id_ordered = [i
                               for i in seq_summary.id_length
                               if i in selected_seq_id]
    curent_id = selected_seq_id_ordered.pop(0)
    with open(args.output.name, 'w') as output_seq_file:
        for seqs in SeqIO.parse(args.input.name, "fasta"):
            if seqs.id == curent_id:
                SeqIO.write(seqs, output_seq_file, "fasta")
                try:
                    curent_id = selected_seq_id_ordered.pop(0)
                except IndexError:
                    break


def main():
    '''
    Create sample of longn reads, sample size is defined as total length
    '''
    args = get_args()
    seq_summary = get_sequences_summary(args.input)
    # check that total length if larger than required sampling
    if seq_summary.total_length < args.total_length:
        print("Input sequences total length: ",
              seq_summary.total_length,
              file=sys.stderr)
        print("Required output total length: ",
              args.total_length,
              file=sys.stderr)
        print("Required sequence length is larger than input data, exiting",
              file=sys.stderr)
        sys.exit(1)

    print("input sequence info")
    print("total length: ", seq_summary.total_length)
    print("number of sequences: ", seq_summary.number_of_sequence)
    make_sequence_sample(args, seq_summary)


if __name__ == "__main__":
    main()
