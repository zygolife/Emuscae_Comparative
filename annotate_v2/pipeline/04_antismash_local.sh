#!/bin/bash
#SBATCH --nodes 1 --ntasks 24 --mem 96G --out logs/antismash.%A.log -J antismash

module load antismash/4.1.0
source activate antismash
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
    CPU=2
fi

DIR=antismash_out
TOPDIR=funannot
FILE=$(ls $DIR/*final.gbk)
INFILE=$(ls $TOPDIR/predict_results/*.gbk)
# removed --asf --cassis
if [ ! $INFILE ]; then
     echo "no file in $TOPDIR/predict_results"
     exit
fi
if [ ! $FILE ]; then
    antismash --taxon fungi --outputfolder $DIR --full-hmmer \
        --clusterblast --smcogs --subclusterblast --knownclusterblast -c $CPU $INFILE
fi
