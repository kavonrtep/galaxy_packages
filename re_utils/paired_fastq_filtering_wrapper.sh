#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# run filtering
WD="`dirname $0`"
ORIDIR=$PWD
cd $WD
WD=$PWD  # absolute path to this script
cd $ORIDIR
SAMPLING=""
TRIM_END=""
TRIM_START=""
PERCENT_ABOVE="95"
CUTADAPT=""
RENAME=""
FILTER_SEQ=""
while getopts "a:b:o:n:c:p:e:s:N:C:G:F:R" OPTION
do
    case $OPTION in
        a)
            FASTAA=$OPTARG;;
        b)
            FASTAB=$OPTARG;;
        o)
            PAIRED_OUTPUT=$OPTARG;;
        n)
            SAMPLING=( -n ${OPTARG} );;
        c)
            CUT_OFF=$OPTARG;;
        p)
            PERCENT_ABOVE=$OPTARG;;
        e)
            TRIM_END=( -e ${OPTARG} );;
        s)
            TRIM_START=( -s ${OPTARG} );;
        N)
            MAX_N=${OPTARG};;
        C)
            CUTADAPT=(-C " "${OPTARG}" " );;
        G)
            PNG_OUTPUT=${OPTARG};;
        R)
            RENAME="-R";;
        F)
            FILTER_SEQ=( -F ${OPTARG} );;

    esac
done
fasta_tmp_fileX=$(mktemp)
fasta_tmp_fileY=$(mktemp)

if [ -z "$CUTADAPT" ] # test if$CUTADAPT is empty
then
    Rscript ${WD}/paired_fastq_filtering.R -a $FASTAA -b $FASTAB -x $fasta_tmp_fileX -y $fasta_tmp_fileY  ${SAMPLING[@]} -c $CUT_OFF -G $PNG_OUTPUT\
         -p $PERCENT_ABOVE  ${TRIM_START[@]}  ${TRIM_END[@]} -N $MAX_N $RENAME ${FILTER_SEQ[@]}
else
    Rscript ${WD}/paired_fastq_filtering.R -a $FASTAA -b $FASTAB -x $fasta_tmp_fileX -y $fasta_tmp_fileY  ${SAMPLING[@]} -c $CUT_OFF -G $PNG_OUTPUT\
         -p $PERCENT_ABOVE  ${TRIM_START[@]}  ${TRIM_END[@]} -N $MAX_N "${CUTADAPT[@]}" $RENAME ${FILTER_SEQ[@]}
fi

python ${WD}/fasta_interlacer.py -a $fasta_tmp_fileX -b $fasta_tmp_fileY -p $PAIRED_OUTPUT -x fasta_tmp_single


rm $fasta_tmp_fileX
rm $fasta_tmp_fileY
