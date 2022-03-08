#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 8gb 

module load hmmer/3

#phmmer --domtbl TcMar1.ORF1.domtbl -E 1e-2 TcMar1.ORF1.aa TcMar/TcMar1.hits.aa.fasta > TcMar1.ORF1.phmmer
#phmmer --domtbl TcMar1.ORF2.domtbl -E 1e-2 TcMar1.ORF2.aa TcMar/TcMar1.hits.aa.fasta > TcMar1.ORF2.phmmer
hmmsearch --cpu 8 --domtbl TcMar1.domtbl -E 1e-3  Mariner.hmm TcMar/TcMar1.hits.aa.fasta > TcMar1.hmmsearch
esl-sfetch --index TcMar/TcMar1.hits.aa.fasta
grep -v '^#' TcMar1.domtbl | awk '{print $1}' | esl-sfetch -f TcMar/TcMar1.hits.aa.fasta - > TcMar1.mariner_hits.aa.fa
hmmalign -o TcMar1.mariner_hits.aa.aln Mariner.hmm TcMar1.mariner_hits.aa.fa  
esl-reformat --replace=x:- --gapsym=- afa TcMar1.mariner_hits.aa.aln | perl -p -e 'if (! /^>/) { s/[ZBzbXx\*]/-/g }' > TcMar1.mariner_hits.aa.clnaln
module load clipkit
clipkit --log TcMar1.mariner_hits.aa.clnaln 
