#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --output=annot_01.%A.log
#SBATCH --time=7-0:00:00
#SBATCH -p batch

module load funannotate/git-live

if [ ! -f config.txt ]; then
 echo "need a config file for parameters"
 exit
fi

source config.txt
if [ ! $GENOME ]; then
 echo "need a genome file defined in config file"
 exit
fi

B=$(basename $GENOME .fasta)
BNAME=$(basename $GENOME)
if [ ! $MIN ]; then
 MIN=500
fi
if [ ! -f $B.clean.fasta ]; then

# mkdir -p /scratch/$USER/$B
# rsync -a -L $GENOME /scratch/$USER/$B
# pushd /scratch/$USER/$B
 funannotate clean -i $BNAME -o $B.clean.fasta -m $MIN
# popd
# rsync -a /scratch/$USER/$B/$B.clean.fasta ./
fi
if [ ! -f $SORTED ]; then
 funannotate sort -i $B.clean.fasta -o $SORTED
fi
