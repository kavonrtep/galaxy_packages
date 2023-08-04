#!/bin/bash

#set -euo pipefail # this is something like string in perl
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
FILTER_SEQ=""
while getopts "a:o:n:c:p:e:s:N:C:G:F:" OPTION
do
    case $OPTION in
        a)
            FASTAA=$OPTARG;;
        o)
            OUTPUT=$OPTARG;;
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
        F)
            FILTER_SEQ=( -F ${OPTARG} );;
    esac
done



if [ -z "$CUTADAPT" ] # test if $CUTADAPT is empty
then
    Rscript ${WD}/single_fastq_filtering.R -a $FASTAA -x $OUTPUT  ${SAMPLING[@]} -c $CUT_OFF\
         -p $PERCENT_ABOVE  ${TRIM_START[@]}  ${TRIM_END[@]} -N $MAX_N -G $PNG_OUTPUT ${FILTER_SEQ[@]}
else
    Rscript ${WD}/single_fastq_filtering.R -a $FASTAA -x $OUTPUT  ${SAMPLING[@]} -c $CUT_OFF -G $PNG_OUTPUT\
         -p $PERCENT_ABOVE  ${TRIM_START[@]}  ${TRIM_END[@]} -N $MAX_N "${CUTADAPT[@]}" ${FILTER_SEQ[@]}
fi




