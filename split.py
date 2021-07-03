#!/usr/bin/env python3

import os
import argparse

parser = argparse.ArgumentParser(description='Split contact pairs into sub folders by chromosome')

parser.add_argument('input',  type=str, help='input contact pairs')
parser.add_argument('--name', type=str, help='folder name', default='.')

args = parser.parse_args()

if args.name != '.':
    if not os.path.exists(args.name):
        os.makedirs(args.name)
if not os.path.exists(args.name+'/intrachr'):
    os.makedirs(args.name+'/intrachr')
if not os.path.exists(args.name+'/interchr'):
    os.makedirs(args.name+'/interchr')


def getChr(ch):
    ch = ch.lstrip('chr')
    if ch[0] == 'X': return 23
    if ch[0] == 'Y': return 24
    if ch[0] == 'M': return 25
    try:
        ch = int(ch)
    except:
        ch = 0
    return ch


outs = {}
with open(args.input, 'r') as f:
    for line in f:
        line = line.rstrip('\n')
        cols = line.split()
        chr1 = cols[1]
        chr2 = cols[6]
        frag1 = cols[5]
        frag2 = cols[10]

        if chr1.find('_') != -1 or chr2.find('_') != -1: continue

        if chr1 == chr2:
            id = 'chr' + chr1.lstrip('chr')
            type = '/intrachr/'
        else:
            c1 = getChr(chr1)
            c2 = getChr(chr2)
            if c1 > c2:
                tmp  = chr1
                chr1 = chr2
                chr2 = tmp
                tmp  = frag1
                frag1 = frag2
                frag2 = tmp
            id = 'chr' + chr1.lstrip('chr') + '_' + chr2.lstrip('chr')
            type = '/interchr/'
        if id not in outs:
            if not os.path.exists(args.name+type+id):
                os.makedirs(args.name+type+id)
            outs[id] = open(args.name+type+id+'/'+id+'.tsv', 'a')
        
        outs[id].write('\t'.join([chr1, frag1, chr2, frag2])+'\n')

for id in outs:
    outs[id].close()



