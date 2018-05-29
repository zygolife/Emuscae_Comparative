#!/usr/bin/bash
#SBATCH --ntasks 2 --nodes 1 --mem 8G --time 1-0:0:0

module load braker/2.0.1
module load augustus/3.2.2
export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.2.2/config
