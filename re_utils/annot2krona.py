#!/usr/bin/env python3
'''
take various inputs and convert it to krona tabular format for visualization
supported inputs:
- DANTE gff3
- TODO PROFREP gff3
- TODO RE archive - normal run
- TODO  RE archive - comparative
-
'''
import argparse
import re
import collections


def parse_dante_gff(f):
    '''load gff3 file and return classification with counts'''
    r = re.compile("Final_Classification=")
    cls_count = collections.defaultdict(int)
    for line in f:
        if re.match("#", line.strip()):
            continue
        attributes = line.split("\t")[8].split(";")
        cls_raw = list(filter(r.match, attributes))[0]
        cls = re.sub(r, "",cls_raw)
        cls_count[cls] += 1

    return cls_count


def export_classification(cls, f):
    '''save classification to tab delimited file'''
    for i in cls:
        f.write('{}\t{}\n'.format(cls[i], i.replace("|","\t")))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--format', choices=['dante', 'profrep', 're'])
    parser.add_argument('-i', '--input', type=argparse.FileType('r'))
    parser.add_argument('-o', '--output', type=argparse.FileType('w'))

    args = parser.parse_args()

    if args.format == "dante":
        classification = parse_dante_gff(args.input)

    if args.format in ["profrep" 're']:
        print("Not implemented")
        exit(0)

    export_classification(classification, args.output)
