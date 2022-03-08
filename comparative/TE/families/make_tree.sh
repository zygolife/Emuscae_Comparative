#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 96 -C xeon --mem 64gb 

module load fasttree

FastTreeMP -wag -gamma < TcMar1.mariner_hits.aa.clnaln.clipkit > TcMar1.mariner_hits.aa.clnaln.clipkit.ft.tre
