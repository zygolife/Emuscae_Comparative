#!/usr/bin/env python3
from Bio import SeqIO
import sys, re

if len(sys.argv) < 2:
    print("need a filename in fasta format as first argument")
    sys.exit()

infile = sys.argv[1]

longestseqs = {}
with open(infile,"r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        genename = re.sub(r'\-T\d+$',"",record.id)
        if genename not in longestseqs or len(record) > len(longestseqs[genename]):
            longestseqs[genename] = record
SeqIO.write(longestseqs.values(), "{}.longest".format(infile), "fasta")
