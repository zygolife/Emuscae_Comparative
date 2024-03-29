#!/bin/bash -l
#SBATCH --nodes 1 -c 24 -n 1 --mem 32G --out logs/busco_euk.%a.log -J busco_euk
module load busco
# for augustus training
#export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.3/config
# set to a local dir to avoid permission issues and pollution in global
export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.3/config)

module load workspace/scratch

CPU=${SLURM_CPUS_ON_NODE}
N=${SLURM_ARRAY_TASK_ID}
if [ ! $CPU ]; then
     CPU=2
fi

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
if [ -z ${SLURM_ARRAY_JOB_ID} ]; then
	SLURM_ARRAY_JOB_ID=$$
fi
GENOMEFOLDER=genomes
LINEAGE=eukaryota_odb10
OUTFOLDER=BUSCO_genome_euk
mkdir -p $OUTFOLDER
SAMPLEFILE=samples.csv
SEED_SPECIES=entomophthora_muscae_ucb
IFS=,
cat $SAMPLEFILE | sed -n ${N}p | while read FILE
do
	BASE=$(basename $FILE .fasta)
	GENOMEFILE=$(realpath $GENOMEFOLDER/$FILE)
	
	if [ -d "$OUTFOLDER/${BASE}" ];  then
	    echo "Already have run $BASE in folder busco - do you need to delete it to rerun?"
	    exit
	else
	    busco -m genome -l $LINEAGE -c $CPU -o ${BASE} --out_path ${OUTFOLDER} --offline --augustus_species $SEED_SPECIES --in $GENOMEFILE --download_path $BUSCO_LINEAGES --tar
	fi
done
