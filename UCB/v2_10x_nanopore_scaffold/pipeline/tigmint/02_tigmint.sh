#!/usr/bin/bash -l
#SBATCH -p intel -N 1 -n 32 --mem 96gb --out logs/tigmint.log
#

# these resources assume the aln files were already made since that runs more 
# efficiently separated among diff nodes

module load tigmint

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi
CPU=$CPUS

mkdir -p tigmint
cd tigmint

ln -s ../genome/emus_flye_assembly_medaka_polish.fasta emus_flye_assembly_medaka_polish.fa
ln -s ../EmusUCB_10X/outs/barcoded.fq.gz EmusUCB_10X.fq.gz
ln -s ../makeTSVfile.py .

tigmint-make tigmint draft=emus_flye_assembly_medaka_polish reads=EmusUCB_10X t=$CPU G=1G
