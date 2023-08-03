#!/usr/bin/env python
import sys
import re
from collections import defaultdict
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-i" ,"--input", type=argparse.FileType('r'), help="path to file CLUSTER_table.csv")
parser.add_argument("-o" ,"--output", type=argparse.FileType('w'), help="output file name")
parser.add_argument("-m", "--use_manual", action='store_true', default=False)

args = parser.parse_args()

column = 6 if args.use_manual else 4
if args.use_manual:
    annotation="Final_annotation"
else:
    annotation="Automatic_annotation"

header = False
clust_info = {}
counts = defaultdict(lambda: 0)
top_clusters = 0
with open(args.input.name, 'r') as f:
    csv_reader = csv.reader(f, delimiter = "\t")
    for parts in csv_reader:
        if len(parts) == 0:
            continue
        if parts[0] == "Cluster" and parts[1]== "Supercluster":
            header = True
            header_columns = parts
            column = header_columns.index(annotation)
            continue
        if header:
            classification = "Top_clusters\t" + "\t".join(parts[column].split("/")[1:]).replace('"','')
            counts[classification] += int(parts[3])
            top_clusters += int(parts[3])
        elif len(parts) >= 2:
            try:
                clust_info[parts[0].replace('"', '')] = int(parts[1])
            except ValueError:
                pass


counts['Singlets'] = clust_info['Number_of_singlets']
counts['Small_cluster'] = int(clust_info['Number_of_reads_in_clusters']) - top_clusters

with open(args.output.name, 'w') as fout:
    for cls_line, nreads in counts.items():
        fout.write(str(nreads) +"\t" + cls_line + "\n")



