#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 32 -C xeon --mem 64gb --out logs/transdecoder.log 

module load transdecoder
module load hmmer/3.3.2-mpi
module load db-pfam

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
    ASM=Emus-Ref3.fasta
    TD=$ASM.transdecoder_dir
    if [ ! -f $TD/longest_orfs.pep ]; then
    	TransDecoder.LongOrfs -t $ASM
    fi
    if [ ! -f $TD/longest_orfs.pep ]; then
	echo "Failed to run transdecoder"
	exit
    fi
    if [ ! -f trinity_$STRAIN.Trinity.hmmscan ]; then
    	srun hmmsearch --mpi --cut_ga --domtbl trinity_$STRAIN.Trinity.domtbl -o trinity_$STRAIN.Trinity.hmmscan $PFAM_DB/Pfam-A.hmm $TD/longest_orfs.pep
    fi
    TransDecoder.Predict -t $ASM --retain_pfam_hits trinity_$STRAIN.Trinity.domtbl
