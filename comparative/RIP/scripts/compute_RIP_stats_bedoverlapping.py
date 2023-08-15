#!/usr/bin/env python

import csv, sys, os, argparse

# genome fai file for sizes
# BED file for RIP locations


parser = argparse.ArgumentParser(description='Compute RIP window summary')
parser.add_argument('--fai',help='The samtools fai index file',required=True,
                    type=argparse.FileType('r'))
parser.add_argument('--bed', type=argparse.FileType('r'),
                    required=True,
                    help='BED file from RIP_index_calculation.pl')

parser.add_argument('--name', required=False,
                    help='Species name')
parser.add_argument('--showheader',required=False,default=False,help='Print header column',action='store_true')
parser.add_argument('--minlen',required=False,default=5000,help='Min length contig to use in calculations')

args = parser.parse_args()

genomesize = 0
RIPwindows = {} # this will be sets initialized assuming RIP positions will be sparse

faiparse = csv.reader(args.fai, delimiter="\t")
for ctg in faiparse:
    ctglen = int(ctg[1])
    if ctglen >= args.minlen:
        genomesize += ctglen
bedparse = csv.reader(args.bed,delimiter="\t")
RIP_length = 0
for bed in bedparse:
    ctg = bed[0]
    start = int(bed[1])-1
    end   = int(bed[2])
    score = float(bed[4])
    if score <= 0:
        continue
    if ctg not in RIPwindows:
        RIPwindows[ctg] = {start}
    for i in range(int(bed[1]),int(bed[2])):
        RIPwindows[ctg].add(i)
for ctg in RIPwindows:
    RIP_length += len(RIPwindows[ctg])
header = ['RIP','Genome_Size','RIP_Percent']
dat    = [RIP_length,genomesize, "%.3f"%(100 * (RIP_length/genomesize))]
if args.name:
    header.append('Name')
    dat.append(args.name)

outcsv = csv.writer(sys.stdout,delimiter="\t")
if args.showheader:
    outcsv.writerow(header)
outcsv.writerow(dat)
