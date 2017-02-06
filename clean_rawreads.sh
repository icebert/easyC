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
    
    if [ $PHRED -eq 33 ]
    then
        FQ_VERSION="-phred33"
    else
        FQ_VERSION="-phred64"
    fi
    
    #
    # When fastqc results are available, try to determine the 
    # phred of fastq through fastqc report.
    #
    ID=${reads1##*/}
    ID=${ID%.gz}
    ID=${ID%.fastq}
    
    ENCODING=""
    if [ -f fastqc/raw/${ID}_fastqc.html ]
    then
        ENCODING=$(cat fastqc/raw/${ID}_fastqc.html | grep "Encoding" \
                   | sed 's/.*<td>Encoding<\/td>//' | cut -f2 -d '>' \
                   | cut -f1 -d'<')
    elif [ -f fastqc/raw/control/${ID}_fastqc.html ]
    then
        ENCODING=$(cat fastqc/raw/control/${ID}_fastqc.html | grep "Encoding" \
                   | sed 's/.*<td>Encoding<\/td>//' | cut -f2 -d '>' \
                   | cut -f1 -d'<')
    fi
    if [ "$ENCODING" == "Sanger / Illumina 1.9" ]
    then
        FQ_VERSION="-phred33"
        if [ $SEPE == "SE" ]
        then
            echo "Detect reads $reads1 base quality are Phred33"
        else
            echo "Detect reads $reads1 and $reads2 base quality are Phred33"
        fi
    elif [ "$ENCODING" == "Illumina 1.3" ]
    then
        FQ_VERSION="-phred64"
        if [ $SEPE == "SE" ]
        then
            echo "Detect reads $reads1 base quality are Phred64"
        else
            echo "Detect reads $reads1 and $reads2 base quality are Phred64"
        fi
    elif [ "$ENCODING" == "Illumina 1.5" ]
    then
        FQ_VERSION="-phred64"
        if [ $SEPE == "SE" ]
        then
            echo "Detect reads $reads1 base quality are Phred64"
        else
            echo "Detect reads $reads1 and $reads2 base quality are Phred64"
        fi
    else
        if [ $SEPE == "SE" ]
        then
            echo "User assume reads $reads1 base quality are Phred$PHRED"
        else
            echo "User assume reads $reads1 and $reads2 base quality are Phred$PHRED"
        fi
    fi
    
    
    if [ $SEPE == "PE" ]
    then
        java -jar $TRIM_HOME/trimmomatic.jar PE $FQ_VERSION -threads $THREAD_NUM \
        "${reads1}" "${reads2}" "${ID1}.clean.fq" "${ID1}.clean.fq.single" \
        "${ID2}.clean.fq" "${ID2}.clean.fq.single" \
        ILLUMINACLIP:$TRIM_HOME/adapters/$ADAPTOR-$SEPE.fa:2:30:7 LEADING:3 \
        TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 1>&2

        if [ $? != 0 ]
        then
            echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tOops, error occurred when \
            clean raw reads\t0\t[FAIL]"
            exit 1
        fi

        rm -f ${ID1}.clean.fq.single ${ID2}.clean.fq.single
    else
        java -jar $TRIM_HOME/trimmomatic.jar SE $FQ_VERSION -threads $THREAD_NUM \
        "${reads1}" "${ID1}.clean.fq" \
        ILLUMINACLIP:$TRIM_HOME/adapters/$ADAPTOR-$SEPE.fa:2:30:7 LEADING:3 \
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



