#!/usr/bin/env python
### this version does not use st input!!!!
# how to use:
# renameSequences.py fasta.file true index.out prefix_length   # for paired sequences
# renameSequences.py fasta.file false index.out prefix_length   # not paired sequences

import sys
paired = sys.argv[2] == "true"
index = open(sys.argv[3], 'w')

if len(sys.argv) == 4:
    prefix = 0
else:
    prefix = int(sys.argv[4])

if paired:
    P = 2
    suffix = "f\n"
else:
    P = 1
    suffix = "\n"
i = j = 0
reader = open(sys.argv[1], mode='r')
for oneline in reader:
    if oneline == "":
        continue
    if oneline[0] == ">":
        i += 1
        j += 1
        prefix_string = oneline[1:(1 + prefix)].strip()
        if j == 1:
            header = ">" + prefix_string + str(i) + suffix
            index.write(oneline[1:].strip() + "\t" + prefix_string + str(i) +
                        suffix)
        if j == 2:
            i -= 1
            header = ">" + prefix_string + str(i) + "r\n"
            index.write(oneline[1:].strip() + "\t" + prefix_string + str(i) +
                        "r\n")
        sys.stdout.write(header)
        if j == P:
            j = 0
    else:
        sys.stdout.write(oneline)
index.close()
