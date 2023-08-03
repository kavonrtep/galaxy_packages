#!/usr/bin/env python
''' interlacing two fastq sequences'''
import sys

def readSingleSeq(file):
    ''' read single seq from fasta file'''
    line = file.readline()
    if not line:
        return False  # end of file
    if line[0] != ">":
        raise Exception("no header on the first line")
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
    ''' write single sequence to fasta file'''
    fileobject.write(">")
    fileobject.write(seq['name'] + "\n")
    fileobject.write(seq['sequence'] + "\n")


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-a",
                      "--fasta_file_A",
                      dest="seqfileA",
                      help="input sequences in fasta format")
    parser.add_option("-b",
                      "--fasta_file_B",
                      dest="seqfileB",
                      help="input sequences in fasta format")
    parser.add_option("-p",
                      "--fasta_file_pairs",
                      dest="seqfile_pairs",
                      help="output file with paired sequences")
    parser.add_option("-x",
                      "--fasta_file_singles",
                      dest="seqfile_singles",
                      help="output file with single sequences")
    options, args = parser.parse_args()

    # Input files
    fA = open(options.seqfileA, 'r')
    fB = open(options.seqfileB, 'r')
    # Output files
    if options.seqfile_pairs:
        fPairs = open(options.seqfile_pairs, 'w')
    else:
        fPairs = open(options.seqfileA + ".pairs", 'w')

    if options.seqfile_singles:
        single = open(options.seqfile_singles, "w")
    else:
        single = open(options.seqfileA + ".single", "w")


    sA1 = readSingleSeq(fA)
    sB1 = readSingleSeq(fB)
    if not sA1 or not sB1:
        raise Exception("\nEmpty sequence on input, nothing to interlace!\n")
    charA = sA1['name'][-1]
    charB = sB1['name'][-1]
    # validate sequence names
    if charA == charB:
        sys.stderr.write(
            "last character of sequence id must be used for distinguishing pairs!")
        exit(1)
        # check first thousand!
    for i in range(3):
        seqA = readSingleSeq(fA)
        seqB = readSingleSeq(fB)
        if (not seqA) or (not seqB):
            # end of file:
            if i == 0:
                sys.stderr.write("input file is empty")
                exit(1)
            else:
                break
        if seqA['name'][-1] == charA and seqB['name'][-1] == charB:
            continue
        else:
            sys.stderr.write(
                "last character of sequence id must be used for distinguishing pairs!")
            exit(1)

    fA.seek(0)
    fB.seek(0)

    buffA = {}
    buffB = {}
    buffA_names = []
    buffB_names = []

    while True:
        seqA = readSingleSeq(fA)
        seqB = readSingleSeq(fB)
        if not seqA and not seqB:
            break  # end of file

        ## validation and direct checking only if not end of files
        if seqA and seqB:
            #validate:
            if not (seqA['name'][-1] == charA and seqB['name'][-1] == charB):
                sys.stderr.write(
                    "last character of sequence id must be used for distinguishing pairs!")
                exit(1)

                # check if current seqs are pairs
            if seqA['name'][:-1] == seqB['name'][:-1]:
                writeSingleSeq(fPairs, seqA)
                writeSingleSeq(fPairs, seqB)
                continue

            ### compare whith buffers
            ### seqA vs buffB
        if seqA:
            if seqA["name"][:-1] in buffB:
                writeSingleSeq(fPairs, seqA)
                seqtmp = {"name": seqA["name"][:-1] + charB,
                          "sequence": buffB[seqA["name"][:-1]]}
                writeSingleSeq(fPairs, seqtmp)
                # can I empty buffA ???
                for i in buffA_names:
                    seqtmp = {"name": i + charA, "sequence": buffA[i]}
                    writeSingleSeq(single, seqtmp)
                buffA = {}
                buffA_names = []

                j = 0
                for i in buffB_names:
                    seqtmp = {"name": i + charB, "sequence": buffB[i]}
                    del buffB[i]
                    j += 1
                    if i == seqA["name"][:-1]:
                        del buffB_names[0:j]
                        break
                    else:
                        writeSingleSeq(single, seqtmp)
            else:
                buffA[seqA["name"][:-1]] = seqA['sequence']
                buffA_names.append(seqA["name"][:-1])

                ### seqA vs buffB
        if seqB:
            if seqB["name"][:-1] in buffA:
                seqtmp = {"name": seqB["name"][:-1] + charA,
                          "sequence": buffA[seqB["name"][:-1]]}
                writeSingleSeq(fPairs, seqtmp)
                writeSingleSeq(fPairs, seqB)
                # can I empty buffB ???
                for i in buffB_names:
                    seqtmp = {"name": i + charB, "sequence": buffB[i]}
                    writeSingleSeq(single, seqtmp)
                buffB = {}
                buffB_names = []

                j = 0
                for i in buffA_names:
                    seqtmp = {"name": i + charA, "sequence": buffA[i]}
                    del buffA[i]
                    j += 1
                    if i == seqB["name"][:-1]:
                        del buffA_names[0:j]
                        break
                    else:
                        writeSingleSeq(single, seqtmp)

            else:
                buffB[seqB["name"][:-1]] = seqB['sequence']
                buffB_names.append(seqB["name"][:-1])

    fA.close()
    fB.close()
    fPairs.close()

    # write rest of singles:
    for i in buffA:
        seqtmp = {"name": i + charA, "sequence": buffA[i]}
        writeSingleSeq(single, seqtmp)
    for i in buffB:
        seqtmp = {"name": i + charB, "sequence": buffB[i]}
        writeSingleSeq(single, seqtmp)
    single.close()


if __name__ == "__main__":
    main()
