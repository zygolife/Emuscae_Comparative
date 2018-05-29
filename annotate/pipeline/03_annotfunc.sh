#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4 --mem 8G
#SBATCH --mem-per-cpu=8G
#SBATCH --output=annot.funannot_03.%A.log
#SBATCH --time=6-0:00:00
#SBATCH -p batch -J annotfunc
module load funannotate/git-live
module load phobius
CPUS=$SLURM_CPUS_ON_NODE

if [ ! $CPUS ]; then
 CPUS=1
fi

if [ ! -f config.txt ]; then
 echo "need a config file for parameters"
 exit
fi

source config.txt
if [ ! $BUSCO ]; then
 BUSCO=fungi_odb9
fi
if [ ! $ODIR ]; then
 ODIR=$(basename `pwd`).'funannot'
fi
MOREFEATURE=""
if [ $TEMPLATE ]; then
 MOREFEATURE="--sbt $TEMPLATE"
fi
funannotate annotate --busco_db $BUSCO -i $ODIR --species "$SPECIES" --strain "$ISOLATE" --cpus $CPUS $EXTRAANNOT $MOREFEATURE
#~/src/funannotate/funannotate annotate --busco_db $BUSCO -i $ODIR --species "$SPECIES" --strain "$ISOLATE" --cpus $CPUS $EXTRAANNOT $MOREFEATURE
