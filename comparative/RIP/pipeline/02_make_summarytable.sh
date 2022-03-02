#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 96 --mem 128gb --out logs/RIPwindow.%a.log -a 1

module load samtools
module load bioperl
module load workspace/scratch
which perl

SAMPLES=samples.csv
INDIR=genomes
INDEXEXT=fai
OUTDIR=RIP_windowsx

cat $SAMPLES | while read FILENAME
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
  NAME=$(echo -n $FILENAME | perl -p -e 's/_/ /g; s/[\._]v\d+\S+//; s/\.fa(asta)//')
  ./scripts/compute_RIP_stats.py --fai $INDEX --bed $OUTFILE --name "$NAME"
done > RIP_stats.tsv
