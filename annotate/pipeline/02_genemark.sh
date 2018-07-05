#!/usr/bin/bash
#SBATCH --ntasks 12 --nodes 1 --mem 96G --time 48:00:00 --out genemark.log

module load genemarkESET

CPUS=$SLURM_CPUS_ON_NODE
TMPDIR=genemark_run
OUTFILE=genemark.gtf
if [ ! $CPUS ]; then
 CPUS=2
fi

if [ ! -f config.txt ]; then
 echo "need a config file for parameters"
 exit
fi

source config.txt
if [ ! $MASKED ]; then
 echo "NEED TO EDIT CONFIG FILE TO SPECIFY THE INPUT GENOME AS VARIABLE: MASKED=GENOMEFILEFFASTA"
 exit
fi

if [ ! $GENEMARK ]; then
 echo "need to edit config and provide GENEMARK variable to point to gmes.mod file"
 exit
fi
if [ ! -d $TMPDIR ]; then
 mkdir -p $TMPDIR
 pushd $TMPDIR
 bp_seqretsplit.pl ../$MASKED
 popd
fi

cmdfile=genemark.jobs
if [ ! -f $cmdfile ]; then
 for file in $TMPDIR/*.fa;
 do 
  scaf=$(basename $file .fa)
  GTF="$TMPDIR/$scaf.gtf"
  if [ ! -f $GTF ]; then
    reformat_fasta.pl --up --soft_mask --native --in $file --out $TMPDIR/$scaf.masked
    unlink $file
    echo "gmhmme3 -s $scaf -f gtf -m $GENEMARK -o $GTF $TMPDIR/$scaf.masked"
  fi 
 done > $cmdfile
fi

parallel -j$CPUS -a $cmdfile 
perl -p -e 'if( ! /^#/ ) { my @row = split(/\t/,$_); my $scaf = $row[0]; $row[-1] =~ s/gene_id\s+"([^"]+)"; transcript_id\s+"([^"]+)"/gene_id "$scaf.$1"; transcript_id "$scaf.$2"/; $_ = join("\t",@row)}' $TMPDIR/*.gtf > genemark.gtf
