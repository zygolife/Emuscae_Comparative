#!/usr/bin/bash

#SBATCH --ntasks 1 --nodes 1 --time 3-0:0:0 -p batch --out antismash_remote.%A.log --mem 24G
source config.txt
module load miniconda2
module load funannotate
funannotate remote -i funannot -m antismash -e $EMAIL

