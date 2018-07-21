#!/bin/bash
#SBATCH --nodes 1 --ntasks 24 --mem 64G -p batch --out mask.%A.log

CPUS=$SLURM_CPUS_ON_NODE
if [ ! $CPUS ]; then
    CPUS=2
fi
module load funannotate/git-live
if [ -f config.txt ]; then
    source config.txt
fi
if [ ! $RMLIB ]; then
    echo "need an RMLIB file created"
    exit
fi
BASE=$(basename $SORTED .fasta)
MASK1=$BASE.rmlib.mask.fasta

if [ ! -f $MASK1 ]; then
funannotate mask -l $RMLIB -i $SORTED -o $MASK1 --cpus $CPUS
fi
MASK2=$BASE.RM_fungi.mask.fasta
if [ ! -f $MASK2 ]; then
funannotate mask -s fungi -i $SORTED -o $MASK2 --cpus $CPUS
fi
if [ ! -f $MASKED ]; then
module load wu-blast
nmerge $MASK1 $MASK2 > $MASKED
fi

