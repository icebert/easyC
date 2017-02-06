#!/bin/env python

import re
import argparse

parser = argparse.ArgumentParser(description='Digest reference genome with restriction enzyme')

parser.add_argument('fasta', type=str, help='genome file in fasta format')

parser.add_argument('--cut', type=str, default='', help='Restriction enzyme cut sequence')

args = parser.parse_args()

id = 1

def cut(chr, seq, cutend):
    global id
    ends = [m.start() for m in re.finditer(cutend, seq)]
    ends.append(len(seq))
    start = 0
    for end in ends:
        if start != end:
            print('\t'.join([chr, str(start), str(end), str(id)]))
            id += 1
        start = end



chr = None
sequence = ''
with open(args.fasta, 'r') as f:
    for line in f:
        line = line.rstrip('\r').rstrip('\n')
        if line == '': continue
        if line[0] == '>':
            if chr:
                cut(chr, sequence, args.cut)
            chr = line.split()[0][1:]
            sequence = ''
        else:
            sequence += line
cut(chr, sequence, args.cut)



