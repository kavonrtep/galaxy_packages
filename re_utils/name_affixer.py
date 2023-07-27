#!/usr/bin/env python
import sys

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--fastq", dest="fastq", help="fastq file")
parser.add_option("-p", "--prefix", dest="prefix", help="prefix to be added to names")
parser.add_option("-s", "--suffix", dest="suffix", help="suffix to be added",default='')
parser.add_option("-n", "--nspace", dest="nspace", help="number of spaces to ignore",default='0')
options, args = parser.parse_args()
nspace=int(options.nspace)

f=open(options.fastq,"r")
j=0
for oneline in f:
    if oneline=="":
        continue
    j+=1
    if j==5:
        j=1
    if not oneline:
        break
   
    if (oneline[0]=="@" and j==1) or (oneline[0]=="+" and len(oneline)>2 and j==3):
        header=" ".join(oneline.split()[:1+nspace])
        header_out=header[0]+options.prefix+header[1:]+options.suffix+"\n"
        sys.stdout.write(header_out)
    else:
        sys.stdout.write(oneline)

