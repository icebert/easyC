#!/bin/bash

#
# This script is to map cleaned paired end reads to reference genome
# and generate the contact pairs.
#

if [ $# -ne 9 ]
then
    echo "Usage: find_contacts.sh <thread number> <remove duplicates> <keep intermediate file>"
    echo "           <reference> <digest genome> <cut end> <min mapping quality>"
    echo "           <Read_1.fq> <Reads_2.fq>"
    echo ""
    exit 1
fi

THREAD_NUM=$1
REMV_DUP=$2
KEEP_FILE=$3
REF_GENOME=$4
DIGEST=$5
CUT_END=$6
MINMQ=$7
shift 7

READS1=$1
READS2=$2

LIB_PATH=$(dirname $(readlink -e $0))


# Create log directory if necessary
if [ ! -d log ]
then
    mkdir log
fi


echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tBegain reads mapping ..."



# if index doesn't exist, first build the index
if [ ! -e "$REF_GENOME.bwt" ]
then
    echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tBegin bwa index\t0\t[OK]"
    
    time bwa index -a bwtsw $REF_GENOME
    
    if [ $? != 0 ]
    then
        echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tOops, error occurred when \
        bwa index\t0\t[FAIL]"
        exit 1
    fi
    echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tbwa index finished\
    successfully\t0\t[OK]"
fi



#
# Map reads to reference genome
#
echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tBegin bwa mem\t0\t[OK]"

ID=$(echo "${READS1}" | sed 's/.*\///' | grep '_' | sed 's/_.*//')
if [ "$ID" == "" ]
then
    ID=$(echo "${READS1}" | rev | cut -d'/' -f2 | rev)
fi

time bwa mem -t $THREAD_NUM "$REF_GENOME" "$READS1" > "${ID}_1.sam"
        
if [ $? -ne 0 ]
then
    echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tOops, error occurred when \
    bwa mem\t0\t[FAIL]"
    exit 2
fi

time bwa mem -t $THREAD_NUM "$REF_GENOME" "$READS2" > "${ID}_2.sam"
        
if [ $? -ne 0 ]
then
    echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tOops, error occurred when \
    bwa mem\t0\t[FAIL]"
    exit 2
fi

echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tbwa mem finished successfully\t0\t[OK]"


echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tReads mapping finished"



#
# Generate contact pairs
#
echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tBegin generating contact pairs\t0\t[OK]"

$LIB_PATH/hicpair.py --minmq $MINMQ --cutend $CUT_END --digest $DIGEST ${ID}_1.sam ${ID}_2.sam \
          1>contacts.tsv 2>log/statistics

if [ $? -ne 0 ]
then
    echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tOops, error occurred when \
    generating contact pairs\t0\t[FAIL]"
    exit 3
fi

cat contacts.tsv | sort -k2,2 -k7,7 -k6,6n -k11,11n -k5,5 -k10,10 -k3,3n -k8,8n --parallel=$THREAD_NUM --buffer-size=8G > contacts.sorted.tsv
rm -f contacts.tsv


if [ "$KEEP_FILE" == "N" ] || [ "$KEEP_FILE" == "n" ]
then
    rm -f ${ID}_1.sam ${ID}_2.sam
fi


if [ "$REMV_DUP" == 'Y' ] || [ "$REMV_DUP" == 'y' ]
then
    $LIB_PATH/removedup.py contacts.sorted.tsv 1>contacts.sorted.rmdup.tsv 2>>log/statistics

    if [ $? -ne 0 ]
    then
        echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tOops, error occurred when \
        remove duplicates of contact pairs\t0\t[FAIL]"
        exit 4
    fi

    rm -f contacts.sorted.tsv
fi


echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tGenerating contact pairs finished successfully\t0\t[OK]"




exit 0


