#!/bin/sh
## same as galaxy terst
./paired_fastq_filtering_wrapper.sh -a test_data/ERR215189_1_part.fastq.gz -b test_data/ERR215189_2_part.fastq.gz -o test_data/single_output.fasta -G test_data/single_output.png -c 10 -N 0 -p 95

echo "paired fastq filtering with defaults"
./paired_fastq_filtering_wrapper.sh -b test_data/ERR215189_2_part.fastq.gz -a test_data/ERR215189_1_part.fastq.gz -o tmp/test2.1.fasta -G tmp/test2.1.png -c 10 -N 0

echo "paired fastq filtering with with sampling"
./paired_fastq_filtering_wrapper.sh -b test_data/ERR215189_2_part.fastq.gz -a test_data/ERR215189_1_part.fastq.gz -o tmp/test2.a2.fasta -G tmp/test2.2.png -c 10 -N 0 -n 500

echo "paired fastq filtering with with sampling"
./paired_fastq_filtering_wrapper.sh -b test_data/ERR215189_2_part.fastq.gz -a test_data/ERR215189_1_part.fastq.gz -o tmp/test2.a3.fasta -G tmp/test2.2.png -c 10 -N 0 -n 653

echo "paired fastq filtering with with sampling"
./paired_fastq_filtering_wrapper.sh -b test_data/ERR215189_2_part.fastq.gz -a test_data/ERR215189_1_part.fastq.gz -o tmp/test2.a4.fasta -G tmp/test2.2.png -c 10 -N 0 -n 547

echo "paired fastq filtering with with sampling"
./paired_fastq_filtering_wrapper.sh -b test_data/ERR215189_2_part.fastq.gz -a test_data/ERR215189_1_part.fastq.gz -o tmp/test2.a5.fasta -G tmp/test2.2.png -c 10 -N 0 -n 839

echo "paired fastq filtering with contaminant removing"
./paired_fastq_filtering_wrapper.sh -F tool_data/organele_ref_and_phi-X174.fasta -b test_data/ERR215189_2_part.fastq.gz -a test_data/ERR215189_1_part.fastq.gz -o tmp/test2.3.fasta -G tmp/test2.3.png -c 10 -N 0

echo "paired fastq filtering with contaminant removing + trimming"
./paired_fastq_filtering_wrapper.sh -F tool_data/organele_ref_and_phi-X174.fasta -b test_data/ERR215189_2_part.fastq.gz -a test_data/ERR215189_1_part.fastq.gz -o tmp/test2.3.fasta -G tmp/test2.3.png -c 10 -N 0 -s 10 -e 80


echo "paired fastq filtering with contaminant removing - full run"
./paired_fastq_filtering_wrapper.sh -F tool_data/organele_ref_and_phi-X174.fasta -b test_data/ERR215189_2.fastq.gz -a test_data/ERR215189_1.fastq.gz -o tmp/test2.3.fasta -G tmp/test2.3.png -c 10 -N 0

