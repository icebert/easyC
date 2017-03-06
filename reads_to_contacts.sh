#!/bin/bash

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
   
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# This script is the integration of all protocols from raw reads processing
# to contact matrix generating.
#
# Written by Wang Meng, 2016-06-13.


# Default settings
LIB_PATH=$(dirname $(readlink -e $0)) #absolute path of this file
REF_GENOME=""                         #reference genome file path
TRIM_ADAPTER="Y"                      #whether to trim adapter given raw reads
THREAD_NUM=6                          #number of threads
PHRED=64                              #the default of base quality is in Phred64
REMV_DUP="Y"                          #whether remove duplicates
ADAPTOR="TruSeq2"                     #the default of adaptor sequence is TruSeq2
KEEP_FILE="N"                         #keep intermediate results
CUT_SEQ=""                            #Restriction enzyme recognize sequence
CUT_POINT=""                          #Restriction enzyme cleavage point
MINMQ=30                              #Minimum mapping quality
WD="./"                               #Working directory


source $LIB_PATH/config

export PATH=$FASTQC_HOME:$BWA_HOME:$SAMTOOLS_HOME:$BEDTOOLS_HOME:$PATH


while getopts "r:xt:c:p:ldka:m:w:" optionName
do
    case "$optionName" in
        r) REF_GENOME="$OPTARG";;
        x) TRIM_ADAPTER="N";;
        t) THREAD_NUM=$OPTARG;;
        c) CUT_SEQ="$OPTARG";;
        p) CUT_POINT=$OPTARG;;
        l) PHRED=33;;
        d) REMV_DUP="N";;
        k) KEEP_FILE="Y";;
        a) ADAPTOR="$OPTARG";;
        m) MINMQ="$OPTARG";;
        w) WD="$OPTARG";;
    esac
done
shift $(($OPTIND - 1))
READS=$@


cd $WD

# Create log directory if necessory
if [ ! -d log ]
then
    mkdir log
fi



#
# Check whether existing enzyme digested genome regions
#
DIGEST=${REF_GENOME%.fasta}
DIGEST=${DIGEST%.fa}.$CUT_SEQ.digest.bed
if [ ! -f $DIGEST ]
then
    echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tDigest genome lost\t0\t[FAIL]" \
    | tee -a log/journal
    echo ""
    echo "------------------------------------------------------------"
    echo "Something bad happened. Details are in the log directory."
    echo "See you next time!"
    echo "------------------------------------------------------------"
    echo ""
    exit 1
fi




#
# Check the quality of input reads 
#
if [ ! -d fastqc/raw ]
then
    echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tCheck quality of input reads begin\t0\t[OK]" \
    | tee -a log/journal

    mkdir -p fastqc/raw
    
    fastqc -o fastqc/raw -t $THREAD_NUM $READS 2>/dev/null 1>&2
    
    rm -f fastqc/raw/*.zip

    echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tCheck quality of input reads finished\t0\t[OK]" \
    | tee -a log/journal
fi


#
# Clean the raw reads
#
if [ $TRIM_ADAPTER == 'Y' ]
then
    $LIB_PATH/clean_rawreads.sh $THREAD_NUM $PHRED PE $TRIM_HOME $ADAPTOR \
    $READS 2>>log/trim.log | tee -a log/journal
    
    STATUS=${PIPESTATUS[0]}
    
    if [ $STATUS -eq 1 ]
    then
        echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tClean raw reads failed\t0\t[FAIL]" \
        | tee -a log/journal
        echo ""
        echo "------------------------------------------------------------"
        echo "Something bad happened. Details are in the log directory."
        echo "See you next time!"
        echo "------------------------------------------------------------"
        echo ""
        exit 1
    fi
    
    PHRED=33
    
    if [ $STATUS -eq 0 ]
    then
        TMP=""
        for item in $READS
        do
            ID=${item##*/}
            ID=${ID%%.*}
            TMP="$TMP ${ID}.clean.fq"
        done
        READS=$TMP
    
        #
        # Check the quality of the cleaned reads
        #
        if [ ! -d fastqc/clean ]
        then
            echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tCheck quality of cleaned reads begin\t0\t[OK]" \
            | tee -a log/journal

            mkdir -p fastqc/clean
            
            fastqc -o fastqc/clean -t $THREAD_NUM $READS 2>/dev/null 1>&2
            
            rm -f fastqc/clean/*.zip

            echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tCheck quality of cleaned reads finished\t0\t[OK]" \
            | tee -a log/journal
        fi
    fi
fi




#
# Find contacts
#
set $READS
while [ $# -gt 0 ]
do
    reads1=$1
    reads2=$2
    shift 2

    $LIB_PATH/find_contacts.sh $THREAD_NUM $REMV_DUP $KEEP_FILE $REF_GENOME \
    $DIGEST $CUT_SEQ $CUT_POINT $MINMQ $reads1 $reads2 2>>log/find_contacts.log | tee -a log/journal

    if [ ${PIPESTATUS[0]} -ne 0 ]
    then
        echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tFind contacts failed\t0\t[FAIL]" \
        | tee -a log/journal
        echo ""
        echo "------------------------------------------------------------"
        echo "Something bad happened. Details are in the log directory."
        echo "See you next time!"
        echo "------------------------------------------------------------"
        echo ""
        exit 1
    fi
done


if [ "$KEEP_FILE" == 'N' ]
then
    # remove cleaned reads
    if [ $TRIM_ADAPTER == 'Y' ]
    then
        rm -f $READS
    fi
fi



exit 0




