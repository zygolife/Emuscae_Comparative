#!/bin/bash -l
#SBATCH -p short -C xeon -N 1 -n 32 --mem 24gb 
module load prodigal
module load parallel
CPU=32
parallel -j $CPU prodigal -a {.}.aa -d {.}.cds -o {.}.prodigal -i {} ::: $(ls *.fas)
cat *.aa > TcMar1.hits.aa.fasta

