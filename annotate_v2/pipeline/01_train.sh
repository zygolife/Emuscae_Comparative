#!/bin/bash
#SBATCH --nodes 1 --ntasks 16 --mem 64G -p intel --out train.%A.log -J trainFun
module unload perl
module load perl/5.22.0
module load python/2.7.12
module load funannotate/git-live
module load augustus/3.3
module load lp_solve
module load genemarkHMM
module load salmon
module unload busco
module load busco/2.0
module load trinity-rnaseq
export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.3/config
export PASAHOME=`dirname $(which Launch_PASA_pipeline.pl)`
export TRINITYHOME
CPUS=$SLURM_CPUS_ON_NODE

if [ ! $CPUS ]; then
    CPUS=2
fi

if [ ! -f config.txt ]; then
    echo "need a config file for parameters"
    exit
fi

source config.txt
if [ ! $SORTED ]; then
    echo "NEED TO EDIT CONFIG FILE TO SPECIFY THE INPUT GENOME AS VARIABLE: SORTED=GENOMEFILEFFASTA"
    exit
fi

if [ ! $ODIR ]; then
     ODIR=$(basename `pwd`)."funannot"
fi

funannotate train -i Entomophthora_muscae_UCB.v2.masked.fasta --trinity Entomophthora_muscae.transcriptome.fasta  --species "Entomopthora muscae" --isolate Berkeley --cpus $CPUS \
    -o funannot --max_intronlen 2000 \
    --left RNAseq/ERR1022665_1.fastq.gz --right RNAseq/ERR1022665_2.fastq.gz
