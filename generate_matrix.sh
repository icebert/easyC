#!/bin/bash

#
# This script is to generate contact matrix
#

if [ $# -ne 7 ]
then
    echo "Usage: generate_matrix.sh <thread number> <name> <resolution>"
    echo "       <normalize> <genome size> <digest genome> <contacts.tsv>"
    echo ""
    exit 1
fi

THREAD_NUM=$1
NAME=$2
RESOLUTION=$3
NORM=$4
GENOME_SIZE=$5
DIGEST=$6
shift 6

LIB_PATH=$(dirname $(readlink -e $0))

export OMP_NUM_THREADS=$THREAD_NUM




rm -rf $NAME/$RESOLUTION

for CON in $@
do
    $LIB_PATH/split.py --name $NAME/$RESOLUTION $CON
done

#
# Process intra-chromosome
#
for chr in $(seq 1 24)
do
    if [ $chr -eq 23 ]
    then
        chr='X'
    elif [ $chr -eq 24 ]
    then
        chr='Y'
    fi
    
    INPUT=$NAME/$RESOLUTION/intrachr/chr$chr/chr$chr.tsv
    if [ -f $INPUT ]
    then
        if [ "$NORM" == 'Y' ] || [ "$NORM" == 'y' ]
        then
            NORM_FILE="${INPUT%.tsv}.KRnorm"
        else
            NORM_FILE=""
        fi

        $LIB_PATH/contactMatrix $RESOLUTION $GENOME_SIZE $DIGEST $INPUT $NORM_FILE > ${INPUT%.tsv}.RAWobserved
        
        if [ $? -ne 0 ]
        then
            echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tOops, error occurred when \
            generate contact matrix for chr$chr\t0\t[FAIL]"
            exit 1
        fi
    fi
    rm -f $INPUT
done

#
# Process inter-chromosome
#
for chr1 in $(seq 1 24)
do
    for chr2 in $(seq $(($chr1+1)) 24)
    do
        if [ $chr1 -eq 23 ]
        then
            chr1='X'
        elif [ $chr1 -eq 24 ]
        then
            chr1='Y'
        fi
        if [ $chr2 -eq 23 ]
        then
            chr2='X'
        elif [ $chr2 -eq 24 ]
        then
            chr2='Y'
        fi

        INPUT=$NAME/$RESOLUTION/interchr/chr${chr1}_${chr2}/chr${chr1}_${chr2}.tsv
        if [ -f $INPUT ]
        then
            if [ "$NORM" == 'Y' ] || [ "$NORM" == 'y' ]
            then
                NORM_FILE="${INPUT%/*}/chr$chr1.KRnorm ${INPUT%/*}/chr$chr2.KRnorm"
            else
                NORM_FILE=""
            fi

            $LIB_PATH/contactMatrix $RESOLUTION $GENOME_SIZE $DIGEST $INPUT $NORM_FILE > ${INPUT%.tsv}.RAWobserved
        
            if [ $? -ne 0 ]
            then
                echo -e "$(date '+%Y-%m-%d %H:%M:%S')\tOops, error occurred when \
                generate contact matrix for chr${chr1}_${chr2}\t0\t[FAIL]"
                exit 1
            fi
        fi
        rm -f $INPUT
    done
done


exit 0


