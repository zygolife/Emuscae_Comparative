#!/bin/bash
#SBATCH -p intel --nodes 1 --ntasks 1 --mem 24G --out logs/fix.%A.log
module load python/2.7.12
module load funannotate/git-live
funannotate fix -i funannot/predict_results/Entomophthora_muscae_Berkeley.gbk -t funannot/predict_results/Entomophthora_muscae_Berkeley.tbl
