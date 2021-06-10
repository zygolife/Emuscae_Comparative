#!/usr/bin/bash
#SBATCH -p short -C xeon -N 1 -n 32 --mem 64gb --out logs/arbitr.log

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
	CPU=$SLURM_CPUS_ON_NODE
fi
module load arbitr
OUTDIR=arbitr
mkdir -p $OUTDIR
GENOME=genome/emus_flye_assembly_medaka_polish.fasta
OUTPREFIX=$OUTDIR/Emus.flye_medaka_polish
ALN=aln/Emus_10X.sortloc.bam

arbitr.py -i $GENOME -o $OUTPREFIX -n $CPU $ALN
