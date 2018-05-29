#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 2 --mem 2G --time 1:00:00

module load phobius
dname=$(basename `pwd`)
phobius -short funannot/predict_results/*.proteins.fa > $dname.phobius.out
