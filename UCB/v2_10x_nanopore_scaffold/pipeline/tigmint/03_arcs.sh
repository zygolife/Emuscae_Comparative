#!/usr/bin/bash -l
#SBATCH -p intel -N 1 -n 32 --mem 96gb --out logs/arcs.log
#

# these resources assume the aln files were already made since that runs more 
# efficiently separated among diff nodes

module load tigmint

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi
CPU=$CPUS

cd tigmint

tigmint-make arcs draft=emus_flye_assembly_medaka_polish reads=EmusUCB_10X t=$CPU G=1G
