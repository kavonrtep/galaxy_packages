# RepeatExplorer utilities #

This repository include utilities for preprocessing of NGS data to suitable format for RepeatExplorer  and TAREAN 
analysis. Each tool include also XML file which define tool interface for Galaxy environment

## Available tools ##

### Paired fastq reads filtering and interlacing ###
tool definition file: `paired_fastq_filtering.xml`

This tool is designed to make memory efficient preprocessing of two fastq files. Output of this file can be used as input of RepeatExplorer clustering. Input files can be in GNU zipped archive (.gz extension). Reads are filtered based on the quality, presence of N bases and adapters. Two input fastq files are procesed in parallel. Only complete pair are kept. As the input files are process in chunks, it is required that pair reads are complete and in the same order in both input files. All reads which pass the quality filter fill be writen into output files. If sampling is specified, only sample of sequences will be returned. Cutadapt us run with this options:

```
--anywhere='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
--anywhere='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
--anywhere='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
--anywhere='ATCTCGTATGCCGTCTTCTGCTTG'
--anywhere='CAAGCAGAAGACGGCATACGAGAT'
--anywhere='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC'
--error-rate=0.05
--times=1 --overlap=15 --discard
```

Order of fastq files processing

1. Trimming (optional)
2. Filter by quality
3. Discard single reads, keep complete pairs
4. Cutadapt filtering
5. Discard single reads, keep complete pairs
6. Sampling (optional)
7. Interlacing two fasta files


### single fastq reads filtering ###
tool definition file: `single_fastq_filtering.xml`

This tool is designed to perform preprocessing
of fastq file. Input files can be in GNU zipped archive (.gz extension). Reads
are filtered based on the quality, presence of N bases and adapters. All reads
which pass the quality filter fill be writen into output files. If sampling is
specified, only sample of sequences will be returned. 

### fasta afixer ###
tool definition file: `fasta_affixer.xml`

Tool for appending prefix and suffix to sequences names in fasta formated sequences. This tool is useful
if you want to do comparative analysis with RepeatExplorer and need to
append sample codes to sequence identifiers

### ChIP-Seq-mapper ###


Analysis of NGS sequences from Chromatin Imunoprecipitation. ChiP and Input reads are mapped to contigs obtained from graph based repetitive sequence clustering to enriched repeats. This method was used in (Neumann et al. 2012). for identification of repetitive sequences associated with cetromeric region.

#### Authors ####
Petr Novak, Jiri Macas, Pavel Neumann, Georg Hermanutz

Biology Centre CAS, Czech Republic


#### Installation and dependencies ####

ChIP-Seq-mapper require NCBI blast to be installed, R programming language with installed R2HTML and base64 packages and python3

#### Usage ####

```
ChipSeqRatioAnalysis.py [-h] [-m MAX_CL] [-n NPROC] -c CHIPSEQ -i
                               INPUTSEQ [-o OUTPUT] [-ht HTML] [-t THRESHOLD]
                               -k CONTIGS

optional arguments:
  -h, --help            show this help message and exit
  -m MAX_CL, --max_cl MAX_CL
                        Sets the maximum cluster number. Default = 200
  -n NPROC, --nproc NPROC
                        Sets the number of cpus to be used. Default = all
                        available
  -c CHIPSEQ, --ChipSeq CHIPSEQ
                        Fasta file containing the Chip Sequences
  -i INPUTSEQ, --InputSeq INPUTSEQ
                        Fasta file containing the Input Sequences
  -o OUTPUT, --output OUTPUT
                        Specify a name for the CSV file to which the output
                        will be save to. Default: ChipSeqRatio.csv                      
  -ht HTML, --html HTML                                                                 
                        Specify a name for the html report. Default :                   
                        ChipSeqRatioReport                                              
  -t THRESHOLD, --threshold THRESHOLD                                                   
                        Optional plot filter. Default: mean ration between              
                        Input hits and Chip hits.                                       
  -k CONTIGS, --Contigs CONTIGS                                                         
                        Contig file for blast 
```

####  References ####
[PLoS Genet. Epub 2012 Jun 21. Stretching the rules: monocentric chromosomes with multiple centromere domains. Neumann P, Navrátilová A, Schroeder-Reiter E, Koblížková A, Steinbauerová V, Chocholová E, Novák P, Wanner G, Macas J.](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002777)


## Dependencies ##

R programming environment with installed packages *optparse* and *ShortRead* (Bioconductor)
python3
cutadapt

## License ##

Copyright (c) 2012 Petr Novak (petr@umbr.cas.cz), Jiri Macas and Pavel Neumann,
Laboratory of Molecular Cytogenetics(http://w3lamc.umbr.cas.cz/lamc/)
Institute of Plant Molecular Biology, Biology Centre AS CR, Ceske Budejovice, Czech Republic

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
