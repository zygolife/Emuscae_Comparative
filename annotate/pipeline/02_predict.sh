#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem 24G
#SBATCH --time=3-00:15:00   
#SBATCH --output=annot_02.%A.log
#SBATCH --job-name="Funnannotate"
module unload perl
module unload perl
module load perl/5.24.0
module load python/2.7.14
module load funannotate/git-live
module load diamond
module load lp_solve
module load genemarkHMM
module load EVM/1.1.1-live
module unload augustus
module unload augustus
module load augustus/3.3
module load minimap2

export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.3/config
CPUS=$SLURM_CPUS_ON_NODE

if [ ! $CPUS ]; then
 CPUS=2
fi

if [ ! -f config.txt ]; then
 echo "need a config file for parameters"
 exit
fi

which python
which augustus
source config.txt
if [ ! $MASKED ]; then 
 echo "NEED TO EDIT CONFIG FILE TO SPECIFY THE INPUT GENOME AS VARIABLE: MASKED=GENOMEFILEFFASTA"
 exit
fi

if [ ! $BUSCO ]; then
 echo "NEED TO PROVIDE A BUSCO FOLDER NAME eg. ascomycota_odb9 fungi_odb9 dikarya_odb9 etc"
 exit
elif [ ! -e $BUSCO ]; then
  ln -s /opt/linux/centos/7.x/x86_64/pkgs/funannotate/share/$BUSCO .
fi
if [ ! -f uniprot_sprot.fasta ]; then
 ln -s /opt/linux/centos/7.x/x86_64/pkgs/funannotate/share/uniprot_sprot.fasta .
fi
if [ ! $PROTEINS ]; then
 PROTEINS=uniprot_sprot.fasta
fi

if [ ! "$EXTRA" ]; then
 EXTRA="--ploidy 1"
fi
if [ ! $ODIR ]; then
 ODIR=$(basename `pwd`)."funannot"
fi

# condition genemark?
CMD="funannotate predict -i $MASKED -s \"$SPECIES\" -o $ODIR --isolate \"$ISOLATE\"  --name $PREFIX --busco_db $BUSCO $AUGUSTUSOPTS \
 --AUGUSTUS_CONFIG_PATH $AUGUSTUS_CONFIG_PATH --transcript_evidence $TRANSCRIPTS --protein_evidence $PROTEINS --cpus $CPUS $EXTRA"
echo $CMD
eval $CMD
