#!/bin/env python

import re
import sys
import argparse

parser = argparse.ArgumentParser(description='Generate HiC contact read pairs')

parser.add_argument('sam1', type=str, help='Mapped reads 1 in sam format')
parser.add_argument('sam2', type=str, help='Mapped reads 2 in sam format')

parser.add_argument('--minmq',    type=int, default=1,  help='The minimum mapping quality [Default 1]')
parser.add_argument('--cutseq',   type=str, default='', help='Restriction enzyme cut sequence')
parser.add_argument('--cutpoint', type=int, default=0, help='Restriction enzyme cleavage point')
parser.add_argument('--digest',   type=str, default='', help='Digested genome in bed format')
parser.add_argument('--out',      type=str, default='stdout', help='output file name [Default stdout]')

args = parser.parse_args()

if args.out == 'stdout':
    out = sys.stdout
else:
    out = open(args.out, 'w')


def complement(nt):
    if nt == 'A' or nt == 'a': return 'T'
    if nt == 'T' or nt == 't': return 'A'
    if nt == 'C' or nt == 'c': return 'G'
    if nt == 'G' or nt == 'g': return 'C'
    return 'N'

def rc(seq):
    cseq = [complement(seq[i]) for i in range(0, len(seq))]
    cseq = ''.join(cseq)
    return cseq[::-1]


def parse(read):
    read = read.split()
    flag = int(read[1])
    if flag & 0x10:
        strand = '-'
    else:
        strand = '+'
    return (read[0], int(read[4]), read[2], int(read[3])-1, read[5], read[9], strand)

def getEnd(start, cigar):
    types = re.split('\d+', cigar)[1:]
    sizes = re.split('[A-Z]', cigar)[:-1]
    sizes = [int(size) for size in sizes]
    length = 0
    for i in range(0, len(types)):
        if types[i] in ('M', 'D', 'N', 'X', 'P'):
            length += sizes[i]
    return start + length

def checkSeq(cigar, seq, cutendl, cutendr):
    types = re.split('\d+', cigar)[1:]
    sizes = re.split('[A-Z]', cigar)[:-1]
    sizes = [int(size) for size in sizes]
    maxs  = 0
    maxi  = -1
    for i in range(0, len(types)-1):
        if (types[i] == 'M' and types[i+1] in ('S', 'H')) or \
           (types[i] in ('S', 'H') and types[i+1] == 'M'):
            if sizes[i] + sizes[i+1] > maxs:
                maxs = sizes[i] + sizes[i+1]
                maxi = i
    if maxi > -1:
        leni = 0
        for i in range(0, maxi+1):
            if types[i] in ('M', 'I', 'S'):
                leni += sizes[i]
        seq1 = seq[0:leni]
        seq2 = seq[leni:]
        seq1 = seq1[-(len(cutendl)*2-1):]
        seq2 = seq2[0:(len(cutendr)*2-1)]
        if seq1 == '':
            end1 = seq2.find(cutendr)
            end2 = seq2.find(rc(cutendl))
            if end1 == -1 and end2 == -1: return False
            if (end1 != -1 and cutendl.endswith(seq2[0:end1])) or (end2 != -1 and rc(cutendr).endswith(seq2[0:end2])):
                return True
            else:
                return False
        if seq2 == '':
            start1 = seq1.find(cutendl)
            start2 = seq1.find(rc(cutendr))
            if start1 == -1 and start2 == -1: return False
            if (start1 != -1 and cutendr.startswith(seq1[(start1+len(cutendl)):])) or (start2 != -1 and rc(cutendl).startswith(seq1[(start2+len(cutendr)):])):
                return True
            else:
                return False
        if types[maxi] == 'M':
            start1 = seq1.find(cutendl)
            start2 = seq1.find(rc(cutendr))
            if start1 == -1 and start2 == -1: return False
            if (start1 != -1 and (seq1[(start1+len(cutendl)):]+seq2).startswith(cutendr)) or (start2 != -1 and (seq1[(start2+len(cutendr)):]+seq2).startswith(rc(cutendl))):
                return True
            else:
                return False
        else:
            end1 = seq2.find(cutendr)
            end2 = seq2.find(rc(cutendl))
            if end1 == -1 and end2 == -1: return False
            if (end1 != -1 and (seq1+seq2[0:end1]).endswith(cutendl)) or (end2 != -1 and (seq1+seq2[0:end2]).endswith(rc(cutendr))):
                return True
            else:
                return False
    return False





sam1 = open(args.sam1, 'r')
sam2 = open(args.sam2, 'r')

digest = {}
if args.digest:
    with open(args.digest, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            chr, start, end, id = line.split()
            start = int(start)
            end   = int(end)
            id    = int(id)
            if chr not in digest:
                digest[chr] = []
            digest[chr] += [id] * (end - start)

# Statistics
rawPairs           = 0
lowQuality         = 0
sameFrag           = 0
incompleteDigest   = 0
crossJunction      = 0
ambigous           = 0
finalContacts      = 0
finalIntraContacts = 0
finalInterContacts = 0


curid  = None
read1  = None
read2  = None
group1 = []
group2 = []
while True:
    if read1 and read2:
        group1 = [read1]
        group2 = [read2]
        id1 = read1[0]
        id2 = read2[0]
        if id1 != id2:
            raise Exception('Both sams should be sorted by reads name')
        curid = id1
    while True:
        read1 = sam1.readline()
        if not read1: break
        if read1[0] == '@': continue
        read1 = parse(read1)
        if curid is None: curid = read1[0]
        if read1[0] == curid:
            group1.append(read1)
        else:
            break
    while True:
        read2 = sam2.readline()
        if not read2: break
        if read2[0] == '@': continue
        read2 = parse(read2)
        if read2[0] == curid:
            group2.append(read2)
        else:
            break

    rawPairs += len(group1) * len(group2)

    #Quality filter
    for i in range(0, len(group1)):
        if group1[i][1] < args.minmq:
            group1[i] = None
            lowQuality += 1
    for i in range(0, len(group2)):
        if group2[i][1] < args.minmq:
            group2[i] = None
            lowQuality += 1

    #Cut end filter
    if args.cutseq:
        if args.cutpoint < len(args.cutseq) - args.cutpoint:
            cutendl = args.cutseq[0:len(args.cutseq)-args.cutpoint]
            cutendr = args.cutseq[args.cutpoint:]
        else:
            cutendl = args.cutseq[0:args.cutpoint]
            cutendr = args.cutseq[len(args.cutseq)-args.cutpoint:]

        if len(group1) > 1:
            for i in range(0, len(group1)):
                if group1[i]:
                    if not checkSeq(group1[i][4], group1[i][5], cutendl, cutendr):
                        group1[i] = None
        if len(group2) > 1:
            for i in range(0, len(group2)):
                if group2[i]:
                    if not checkSeq(group2[i][4], group2[i][5], cutendl, cutendr):
                        group2[i] = None

    #Output
    output = []
    for r1 in group1:
        if not r1: continue
        for r2 in group2:
            if not r2: continue
            frag1 = 0
            frag2 = 0
            if digest:
                frag1 = digest[r1[2]][r1[3]]
                frag2 = digest[r2[2]][r2[3]]
                if frag1 == frag2:
                    sameFrag += 1
                    continue
                if len(group1) > 1 or len(group2) > 1:
                    if abs(frag1 - frag2) <= 3:
                        incompleteDigest += 1
                        continue
                    crossJunction += 1
            output.append((r1, r2, frag1, frag2))
    if len(output) == 2:
        (r11, r12, frag11, frag12) = output[0]
        (r21, r22, frag21, frag22) = output[1]
        if frag11 == frag22 and frag12 == frag21:
            del output[1]
    if len(output) == 1:
        (r1, r2, frag1, frag2) = output[0]
        finalContacts += 1
        if r1[2] == r2[2]:
            finalIntraContacts += 1
        else:
            finalInterContacts += 1
        id = [ r1[0] ]
        p1 = [ r1[2], r1[3], getEnd(r1[3], r1[4]), r1[6], frag1 ]
        p2 = [ r2[2], r2[3], getEnd(r2[3], r2[4]), r2[6], frag2 ]
        if (p1[0] > p2[0]) or (p1[0] == p2[0] and p1[1] > p2[1]):
            tmp = p1
            p1 = p2
            p2 = tmp
            if p1[3] == '+':
                p1[3] = '-'
            else:
                p1[3] = '+'
            if p2[3] == '+':
                p2[3] = '-'
            else:
                p2[3] = '+'
        p1 = [str(item) for item in p1]
        p2 = [str(item) for item in p2]
        out.write('\t'.join(id + p1 + p2) + '\n')
    elif len(output) > 1:
        ambigous += len(output)

    if (not read1) or (not read2):
        break

sam1.close()
sam2.close()
if args.out != 'stdout':
    out.close()

#Output statistics
sys.stderr.write('Raw Pairs:\t'                       + str(rawPairs)           + '\n')
sys.stderr.write('Low Quality Pairs:\t'               + str(lowQuality)         + '\n')
sys.stderr.write('Same Fragment Pairs:\t'             + str(sameFrag)           + '\n')
sys.stderr.write('Incomplete Enzyme Digest:\t'        + str(incompleteDigest)   + '\n')
sys.stderr.write('Ambiguous Cross Junction Pairs:\t'  + str(ambigous)           + '\n')
sys.stderr.write('Cross Junction Valid:\t'            + str(crossJunction - ambigous) + '\n')
sys.stderr.write('Final Contacts:\t'                  + str(finalContacts)            + '\n')
sys.stderr.write('Final Contacts Intra-chromosome:\t' + str(finalIntraContacts)       + '\n')
sys.stderr.write('Final Contacts Inter-chromosome:\t' + str(finalInterContacts)       + '\n')







