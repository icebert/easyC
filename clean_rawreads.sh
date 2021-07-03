#!/bin/bash

#
# This script is to clean the raw reads to remove adaptors sequences
# and low quality leading and tailing bases. The input is:
# clean_rawreads.sh <threads> <phred> <single end or paired end> 
# <Trimmomatic home> <adaptor> <file1> ... <fileN>. 
#


if [ $# -lt 6 ]
then
    echo "Usage: clean_rawreads.sh <thread number> <Phred> <SE or PE>"
    echo "           <trimmomatic home> <adaptor> <file1> ... <fileN>"
    echo ""
    exit 1
fi



THREAD_NUM=$1
PHRED=$2
SEPE=$3
TRIM_HOME=$4
ADAPTOR=$5
shift 5

INPUT=$@


echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tClean raw reads begin\t0\t[OK]"


READS=""
set $INPUT
while [ $# -gt 0 ]
do
    item1=$1
    shift
    item2=""
    if [ $SEPE == "PE" ]
    then
        item2=$1
        shift
    fi
    
    ID=${item1##*/}
    ID=${ID%%.*}
    
    if [ -e ${ID}.clean.fq ]
    then
        if [ $SEPE == "PE" ]
        then
            echo "Reads $item1 and $item2 had been cleaned, skip"
        else
            echo "Reads $item1 had been cleaned, skip"
        fi
    else
        READS="$READS $item1 $item2"
    fi
done


READS=$(echo $READS | sed 's/^ *//')
if [ "$READS" == "" ]
then
    echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tNothing to clean, exit cleaning"
    exit 0
fi




set $READS
while [ $# -gt 0 ]
do
    reads1=$1
    shift
    
    ID1=${reads1##*/}
    ID1=${ID1%%.*}
    
    if [ $SEPE == "PE" ]
    then
        reads2=$1
        shift
        
        ID2=${reads2##*/}
        ID2=${ID2%%.*}
    fi
    
    
    if [ $SEPE == "PE" ]
    then
        java -jar $TRIM_HOME/trimmomatic.jar PE -threads $THREAD_NUM \
        "${reads1}" "${reads2}" "${ID1}.clean.fq" "${ID1}.clean.fq.single" \
        "${ID2}.clean.fq" "${ID2}.clean.fq.single" \
        ILLUMINACLIP:$TRIM_HOME/adapters/$ADAPTOR-$SEPE.fa:2:30:10:1:true LEADING:3 \
        TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 1>&2

        if [ $? != 0 ]
        then
            echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tOops, error occurred when \
            clean raw reads\t0\t[FAIL]"
            exit 1
        fi

        rm -f ${ID1}.clean.fq.single ${ID2}.clean.fq.single
    else
        java -jar $TRIM_HOME/trimmomatic.jar SE -threads $THREAD_NUM \
        "${reads1}" "${ID1}.clean.fq" \
        ILLUMINACLIP:$TRIM_HOME/adapters/$ADAPTOR-$SEPE.fa:2:30:10 LEADING:3 \
        TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 1>&2

        if [ $? != 0 ]
        then
            echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tOops, error occurred when \
            clean raw reads\t0\t[FAIL]"
            exit 1
        fi
    fi
done

echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tClean raw reads finished\t0\t[OK]"

exit 0



