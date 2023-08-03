#!/usr/bin/env python3
''' fasta affixer - adding prefixes and suffixes to fasta sequence names'''
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", type=str, help="fasta file")
parser.add_argument("-o", "--output", type=str, help="output fasta file")
parser.add_argument(
    "-p", "--prefix",
    type=str, help="prefix to be added to names")
parser.add_argument(
    "-s", "--suffix",
    type=str, help="suffix to be added",
    default='')
parser.add_argument("-n",
                    "--nspace",
                    type=int,
                    help="number of spaces to ignore",
                    default='0')

args = parser.parse_args()

with open(args.fasta, "r") as f, open(args.output, "w") as out:
    for oneline in f:
        if oneline == "":
            continue
        if not oneline:
            break
        if oneline[0] == ">":
            header = " ".join(oneline.split()[:1 + args.nspace])
            header_out = header[0] + args.prefix + header[1:] + args.suffix + "\n"
            out.write(header_out)
        else:
            out.write(oneline)

