#!/usr/bin/bash -l
#SBATCH --time 5-0:0:0 -p intel -N 1 -n 32 --mem 128gb --out logs/orthofinder_treebuild.%A.log
ulimit -Sn
ulimit -Hn
ulimit -n 80000
ulimit -Sn
ulimit -Hn
CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
	CPU=$SLURM_CPUS_ON_NODE
fi
mkdir -p logs
module load orthofinder/2.5.2
export TMPDIR=/scratch
orthofinder -ft OrthoFinder_diamond/Blast_results/OrthoFinder/Results_Aug19 -t $CPU -a $CPU -S diamond_ultra_sens
