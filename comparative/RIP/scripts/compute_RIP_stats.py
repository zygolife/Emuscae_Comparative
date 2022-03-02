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

args = parser.parse_args()

chromsizes = {}
genomesize = 0

faiparse = csv.reader(args.fai, delimiter="\t")
for ctg in faiparse:
    chromsizes[ctg[0]] = ctg[1]
    genomesize = int(ctg[2]) # last writer wins FAI file 3rd column is accumulated genome size
bedparse = csv.reader(args.bed,delimiter="\t")
RIPwindows = {}
RIP_length = 0
for bed in bedparse:
    ctg = bed[0]
    if ctg not in RIPwindows:
        RIPwindows[ctg] = []
    RIPwindows[ctg].append(bed)
    RIP_length += abs(int(bed[2]) - int(bed[1])+1)

header = ['RIP','Genome_Size','RIP_Percent']
dat    = [RIP_length,genomesize, "%.3f"%(100 * (RIP_length/genomesize))]
if args.name:
    header.append('Name')
    dat.append(args.name)

outcsv = csv.writer(sys.stdout,delimiter="\t")
outcsv.writerow(header)
outcsv.writerow(dat)
