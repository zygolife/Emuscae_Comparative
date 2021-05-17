#!/usr/bin/bash
#SBATCH -p short -N 1 -n 16 --mem 96gb --out logs/map_10x_combine.log
module load bwa
module load samtools/1.11
CPUS=$SLURM_CPUS_ON_NODE
ALN=aln
GENOME1=genome/bwa1/emus_flye_assembly_medaka_polish.fasta
GENOME2=genome/emus_flye_assembly_medaka_polish.fasta
if [ -z $CPUS ]; then
 CPUS=1
fi
CPU=$CPUS
date
hostname

NAME=Emus_10x_ONT
EXT=sort_bx.bam
OUT=$ALN/$NAME.sortbx.bam
# Determine CPUS
if [[ -z ${SLURM_CPUS_ON_NODE} ]]; then
    CPUS=1
else
    CPUS=${SLURM_CPUS_ON_NODE}
fi

if [ ! -s $OUT ]; then
	#which bwa-mem2
	#bwa-mem2 mem -C -t $CPU -p $GENOME2 $READS | samtools view -u -F4 > $TEMP/$NAME.sam
	samtools merge -h $ALN/$NAME.00.$EXT -n --threads $CPU $OUT $ALN/$NAME.*.$EXT
fi
