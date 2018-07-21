#!/bin/bash
#SBATCH --ntasks 48 --nodes 1 --mem 128G -p intel --time 36:00:00 --out iprscan.%A.log
module load python/2.7.12
module load funannotate/git-live
module load iprscan
which python

CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
            CPU=2
fi

FUNOUT=funannot
XML=$FUNOUT/annotate_misc/iprscan.xml
PEP=$FUNOUT/annotate_misc/genome.proteins.fasta
if [ ! -f $XML ]; then
    funannotate iprscan -i $FUNOUT -m local -c $CPU --iprscan_path $(which interproscan.sh)
fi
