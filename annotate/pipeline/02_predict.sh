#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem 24G
#SBATCH --time=3-00:15:00   
#SBATCH --output=annot_02.%A.log
#SBATCH --job-name="Funnannotate"
module unload miniconda2
module load python/2.7.14
module load diamond
module load funannotate/git-live
module load lp_solve
module load genemarkHMM
module unload salmon
module load salmon/0.8.2
module unload busco
module load busco/2.0
module load EVM
module unload augustus
module unload augustus
module load augustus/3.3

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
if [ ! $SORTED ]; then 
 echo "NEED TO EDIT CONFIG FILE TO SPECIFY THE INPUT GENOME AS VARIABLE: SORTED=GENOMEFILEFFASTA"
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
CMD="funannotate predict -i $SORTED -s \"$SPECIES\" -o $ODIR --isolate \"$ISOLATE\"  --name $PREFIX --busco_db $BUSCO --genemark_mod $GENEMARK $AUGUSTUSOPTS \
 --AUGUSTUS_CONFIG_PATH $AUGUSTUS_CONFIG_PATH --transcript_evidence $TRANSCRIPTS --protein_evidence $PROTEINS --cpus $CPUS $EXTRA"
echo $CMD

funannotate predict -i $SORTED -s "$SPECIES" -o $ODIR --isolate "$ISOLATE"  --name $PREFIX --busco_db $BUSCO $AUGUSTUSOPTS --genemark_mod $GENEMARK --AUGUSTUS_CONFIG_PATH $AUGUSTUS_CONFIG_PATH --transcript_evidence $TRANSCRIPTS --protein_evidence $PROTEINS --cpus $CPUS $EXTRA
