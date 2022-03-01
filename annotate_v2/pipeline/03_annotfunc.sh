#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4 --mem 8G
#SBATCH --mem-per-cpu=8G
#SBATCH --output=annot.funannot_03.%A.log
#SBATCH --time=6-0:00:00
#SBATCH -p batch -J annotfunc
module unload miniconda2
module unload miniconda3
module load perl/5.24.0
module load python/2.7.12
module load funannotate/git-live
module switch augustus/3.3
module load phobius

export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.3/config

CPUS=$SLURM_CPUS_ON_NODE

if [ ! $CPUS ]; then
 CPUS=1
fi
EXTRAFEATURE=""
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
if [ $TEMPLATE ]; then
 EXTRAFEATURE="--sbt $TEMPLATE"
fi
funannotate annotate --busco_db $BUSCO -i $ODIR --species "$SPECIES" --strain "$ISOLATE" --cpus $CPUS $EXTRAANNOT --tbl2asn '-l paired-ends' $EXTRAFEATURE
