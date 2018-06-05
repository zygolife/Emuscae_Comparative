#!/bin/bash
#SBATCH --nodes 1 --ntasks 2 --mem 4G --out maker.%A_%a.log 
module unload miniconda2
module unload python
module unload perl
module unload perl
module load perl/5.24.0
module load maker/2.31.8
module unload augustus
module unload augustus
module unload augustus
module load augustus/3.3
which augustus
export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.3/config
maker 
