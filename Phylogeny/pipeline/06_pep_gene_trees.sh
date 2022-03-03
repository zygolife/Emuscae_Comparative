#!/usr/bin/bash
#SBATCH -p short -n 48 --mem 32gb  -N 1 --out logs/make_pep_trees.%A.log
module load fasttree
module load IQ-TREE/2.1.1

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
make -f PHYling_unified/util/makefiles/Makefile.trees HMM=fungi_odb10 -j $CPU PEP

