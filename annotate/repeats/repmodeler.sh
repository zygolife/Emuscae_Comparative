#!/bin/bash
#SBATCH --nodes 1 --ntasks 32 -p intel --mem 64G --time 72:00:00 --out rmodeler.%A.log

module load RepeatModeler
CPU=32
RepeatModeler -pa $CPU -engine ncbi -database Entomophthora_muscae_UCB.500k
