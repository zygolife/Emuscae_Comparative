#!/usr/bin/bash
#SBATCH -p intel,batch -N 1 -n 24 --mem 64gb --out logs/demux_longranger.log

module load longranger
GENOME=genome/emus_flye_assembly_medaka_polish.fasta
CPU=24
MEM=64
#longranger basic --fastqs=10x_reads --localcores=$CPU --localmem=$MEM --mode=local --id=EmusUCB_10X
longranger basic --fastqs=10x_reads --id=EmusUCB_10X  --localmem=$MEM --localcores=$CPU
