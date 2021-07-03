#!/usr/bin/env python3

import re
import argparse

parser = argparse.ArgumentParser(description='Digest reference genome with restriction enzyme')

parser.add_argument('fasta', type=str, help='genome file in fasta format')

parser.add_argument('--seq', type=str, default='', help='Restriction enzyme cut sequence')
parser.add_argument('--point', type=int, default=0, help='Restriction enzyme cleavage point')

args = parser.parse_args()

id = 1


def complement(nt):
    if nt == 'A' or nt == 'a': return 'T'
    if nt == 'T' or nt == 't': return 'A'
    if nt == 'C' or nt == 'c': return 'G'
    if nt == 'G' or nt == 'g': return 'C'
    return 'N'

def reverse_complement(seq):
    cseq = [complement(seq[i]) for i in range(0, len(seq))]
    cseq = ''.join(cseq)
    return cseq[::-1]


def cut(chr, seq, cutseq, cutpoint):
    global id
    rcutseq = reverse_complement(cutseq)
    if rcutseq == cutseq:
        ends = [m.start() + cutpoint for m in re.finditer(cutseq, seq)]
    else:
        ends1 = [m.start() + cutpoint for m in re.finditer(cutseq, seq)]
        ends2 = [m.start() + cutpoint for m in re.finditer(rcutseq, seq)]
        ends  = sorted(set(ends1 + ends2))
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
                cut(chr, sequence, args.seq, args.point)
            chr = line.split()[0][1:]
            sequence = ''
        else:
            sequence += line
cut(chr, sequence, args.seq, args.point)



