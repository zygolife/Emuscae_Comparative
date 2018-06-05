#!/usr/bin/bash 

#SBATCH --nodes=1
#SBATCH --ntasks=64 
#SBATCH  --mem 128G -p stajichlab
#SBATCH  --time=48:00:00
#SBATCH --job-name genemark
#SBATCH --output=Genemarkhmm.%A.out

module load genemarkESET
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
echo "CPU is $CPU"
LINE=${SLURM_ARRAY_TASK_ID}
SPECIES=Entomophthora_muscae_UCB
echo $SPECIES
if [ ! -f $SPECIES.long.fasta ]; then
perl /rhome/jstajich/src/genome-scripts/gene_prediction/select_long_ctgs.pl -m 100000 $SPECIES.fasta > $SPECIES.long.fasta
fi
nohup gmes_petap.pl --min_contig 100000 --cores 64 --fungus --ES --sequence $SPECIES.long.fasta >& train.log
