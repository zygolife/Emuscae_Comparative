#!/bin/bash
#SBATCH --nodes 1 --ntasks 8 --mem 16G --out logs/antismash.%a.log -J antismash

module unload miniconda2
module unload miniconda3
module unload anaconda3
module load antismash/6.0.0
module load antismash/6.0.0
which perl
which antismash
hostname
CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
GBK=$(ls *.gbff)
    #	antismash --taxon fungi --output-dir $OUTDIR/$name/antismash_local  --genefinding-tool none \
      #    --asf --fullhmmer --cassis --clusterhmmer --asf --cb-general --pfam2go --cb-subclusters --cb-knownclusters -c $CPU \
      #    $OUTDIR/$name/$INPUTFOLDER/*.gbk
    time antismash --taxon fungi --output-dir antismash_local \
      --genefinding-tool none --fullhmmer --clusterhmmer --cb-general  --cassis --asf --cb-subclusters --cb-knownclusters \
      --pfam2go -c $CPU $GBK
