#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 96 --mem 128gb --out logs/RIPwindow.%a.log -a 1

module load samtools
module load bioperl
module load workspace/scratch
which perl

CPUS=${SLURM_CPUS_ON_NODE}
if [ ! $CPUS ]; then
 CPUS=1
fi

N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
 N=$1
 if [ -z $N ]; then
   N=1
   echo "defaulting to N value is 1 - specify with --array or cmdline"
 fi
fi

runRIP() {
  IN=$1
  OUTF=$OUTTMP/${IN}.bed
  if [ ! -f $OUTF ]; then
    samtools faidx $INFILE $IN | perl $SCRIPT -r bed -o $OUTF -
  fi
}
export -f runRIP

SCRIPT=scripts/RIP_index_calculation.pl
export SCRIPT
SAMPLES=samples.csv
INDIR=genomes
INDEXEXT=fai
OUTDIR=RIP_windows
mkdir -p $OUTDIR

sed -n ${N}p $SAMPLES | while read FILENAME
do
  INFILE=$INDIR/$FILENAME
  export INFILE
  INDEX=$INFILE.${INDEXEXT}
  if [ ! -f $INDEX ]; then
      samtools faidx $INFILE
  fi
  BASE=$(echo -n $FILENAME | perl -p -e 's/\.fa(sta)?//; s/\.(masked|sorted)//; ')
  OUTTMP=$OUTDIR/$BASE.d
  export OUTTMP
  mkdir -p $OUTTMP
  awk '$2 > 5000' $INDEX | cut -f1 | parallel -j $CPUS runRIP
done
