#!/usr/bin/env python3

import sys
import argparse

parser = argparse.ArgumentParser(description='Generated bed file for binned genome at given resolution')

parser.add_argument('--gsize', type=str, help='genome size file')
parser.add_argument('--chr',   type=str, help='chromosome')
parser.add_argument('--res',   type=int, help='resolution')

args = parser.parse_args()

ch = args.chr.lstrip('chr')
length = 0
with open(args.gsize, 'r') as f:
    for line in f:
        c, size = line.rstrip('\n').split()
        c = c.lstrip('chr')
        if c == ch:
            length = int(size)
            break
if length == 0:
    sys.stderr.write('Cannot find {0} in genome size file {1}\n\n'.format(args.chr, args.gsize))

start = 0
while start < length:
    sys.stdout.write('\t'.join([args.chr, str(start), str(start + args.res), str(start)]) + '\n')
    start += args.res

