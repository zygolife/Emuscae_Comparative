#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 16 --mem 96G --time 12:00:00 --out emapper.log

module load eggnog-mapper
dname=$(basename `pwd`)
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
 CPU=2
fi
emapper.py -i funannot/predict_results/*.proteins.fa -d fuNOG -o $dname --cpu $CPU --usemem -m diamond 
