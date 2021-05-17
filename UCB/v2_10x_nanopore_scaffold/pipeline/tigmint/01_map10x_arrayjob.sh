#!/usr/bin/bash
#SBATCH -p short -N 1 -n 24 --mem 32gb --out logs/map_10x.%a.log
module load bwa
module load bwa-mem2/2.2.1
module load samtools/1.11
CPUS=$SLURM_CPUS_ON_NODE
ALN=aln
GENOME1=genome/bwa1/emus_flye_assembly_medaka_polish.fasta
GENOME2=genome/emus_flye_assembly_medaka_polish.fasta
mkdir -p $ALN
if [ -z $CPUS ]; then
 CPUS=1
fi
CPU=$CPUS
date
hostname

NAME=Emus_10x_ONT
EXT=fq.gz
FWD=EmusUCB_10X/outs/split/bc_R1_
REV=EmusUCB_10X/outs/split/bc_R2_

TEMP=/scratch/$NAME.$$

# Determine CPUS
if [[ -z ${SLURM_CPUS_ON_NODE} ]]; then
    CPUS=1
else
    CPUS=${SLURM_CPUS_ON_NODE}
fi


N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi

N=$(printf %02d $N)

if [ ! -s $ALN/$NAME.$N.sort_bx.bam ]; then
	mkdir -p $TEMP
	#which bwa-mem2
	#bwa-mem2 mem -C -t $CPU -p $GENOME2 $READS | samtools view -u -F4 > $TEMP/$NAME.sam
	FWDFILE=${FWD}${N}.${EXT}
	REVFILE=${REV}${N}.${EXT}
	bwa-mem2 mem -C -t $CPU $GENOME2 $FWDFILE $REVFILE > $TEMP/$NAME.$N.sam
	samtools sort -T $TEMP/$NAME -m 1G --threads $CPU -tBX -o $ALN/$NAME.$N.sort_bx.bam $TEMP/$NAME.$N.sam
	rm -rf $TEMP
fi
