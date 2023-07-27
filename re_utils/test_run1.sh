#!/bin/sh

# same like galaxy test:

./single_fastq_filtering_wrapper.sh -a test_data/ERR215189_1_part.fastq.gz -o test_data/single_output.fasta -G test_data/single_output.png -c 10 -N 0 -p 95


echo "single fastq filtering with defaults"
./single_fastq_filtering_wrapper.sh -a test_data/ERR215189_1_part.fastq.gz -o tmp/test1.fasta -G tmp/test1.png -c 10 -N 0

echo "single fastq filtering with with sampling"
./single_fastq_filtering_wrapper.sh -a test_data/ERR215189_1_part.fastq.gz -o tmp/test2a.fasta -G tmp/test2a.png -c 10 -N 0 -n 500

echo "single fastq filtering with with sampling"
./single_fastq_filtering_wrapper.sh -a test_data/ERR215189_1_part.fastq.gz -o tmp/test2b.fasta -G tmp/test2b.png -c 10 -N 0 -n 647

echo "single fastq filtering with with sampling"
./single_fastq_filtering_wrapper.sh -a test_data/ERR215189_1_part.fastq.gz -o tmp/test2c.fasta -G tmp/test2c.png -c 10 -N 0 -n 839

echo "single fastq filtering with with sampling"
./single_fastq_filtering_wrapper.sh -a test_data/ERR215189_1_part.fastq.gz -o tmp/test2d.fasta -G tmp/test2d.png -c 10 -N 0 -n 911


echo "single fastq filtering with contaminant removing"
./single_fastq_filtering_wrapper.sh -F tool_data/organele_ref_and_phi-X174.fasta -a test_data/ERR215189_1_part.fastq.gz -o tmp/test3.fasta -G tmp/test3.png -c 10 -N 0 


