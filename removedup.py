#!/usr/bin/env python3

import sys
import argparse

parser = argparse.ArgumentParser(description='Remove duplicates in HiC contact pairs')

parser.add_argument('input',    type=str, help='input contact pairs')
parser.add_argument('--within', type=int, default=4, help='The max distance of two positions being considered as duplicates [Default 4]')
parser.add_argument('--out',    type=str, default='stdout', help='output file name [Default stdout]')

args = parser.parse_args()

if args.out == 'stdout':
    out = sys.stdout
else:
    out = open(args.out, 'w')



def isDup(a, b):
    a1_chr = a[1]
    a1_pos = int(a[2])
    a1_strand = a[4]
    a2_chr = a[6]
    a2_pos = int(a[7])
    a2_strand = a[9]
    b1_chr = b[1]
    b1_pos = int(b[2])
    b1_strand = b[4]
    b2_chr = b[6]
    b2_pos = int(b[7])
    b2_strand = b[9]
    if a1_chr == b1_chr and a1_strand == b1_strand and abs(a1_pos - b1_pos) <= args.within and \
       a2_chr == b2_chr and a2_strand == b2_strand and abs(a2_pos - b2_pos) <= args.within:
        return True
    return False




duplicates = 0

with open(args.input, 'r') as input:
    last = None
    for line in input:
        if not last:
            out.write(line)
            line = line.rstrip('\n')
            last = line.split()
        else:
            nline = line.rstrip('\n')
            cur   = nline.split()
            if isDup(last, cur):
                duplicates += 1
            else:
                out.write(line)
                last = cur

if args.out != 'stdout':
    out.close()

#Output statistics
sys.stderr.write('Dupicates:\t' + str(duplicates) + '\n')



