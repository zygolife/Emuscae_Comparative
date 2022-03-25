#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 2 --mem 16gb --out logs/RIP_summary.log 

module load samtools
module load bioperl
module load workspace/scratch
which perl

SAMPLES=samples.csv
INDIR=genomes
INDEXEXT=fai
OUTDIR=RIP_windows
FIRST=1
cat $SAMPLES | grep -v NM02 | while read FILENAME
do
  INFILE=$INDIR/$FILENAME
  INDEX=$INFILE.${INDEXEXT}
  if [ ! -f $INDEX ]; then
      samtools faidx $INFILE
  fi
  BASE=$(echo -n $FILENAME | perl -p -e 's/\.fa(sta)?//; s/\.(masked|sorted)//; ')
  OUTTMP=$OUTDIR/$BASE.d
  OUTFILE=$OUTDIR/$BASE.RIP.bed
  cat $OUTTMP/*.bed > $OUTFILE
  NAME=$(echo -n $FILENAME | perl -p -e 's/_/ /g; s/[\._]v\d+\S+//; s/\.fa(sta)//')
  if [ ! -z $FIRST ]; then
    ./scripts/compute_RIP_stats_bedoverlapping.py --fai $INDEX --bed $OUTFILE --name "$NAME" --showheader
  else
    ./scripts/compute_RIP_stats_bedoverlapping.py --fai $INDEX --bed $OUTFILE --name "$NAME"
  fi
  FIRST=""
done > RIP_stats.tsv
