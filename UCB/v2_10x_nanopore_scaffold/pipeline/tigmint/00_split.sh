#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 4 --mem 16gb --out logs/split_fastq.log

#module load fastqsplitter

mkdir -p EmusUCB_10X/outs/split
pushd EmusUCB_10X/outs/split

if [[ -z ${SLURM_CPUS_ON_NODE} ]]; then
    CPUS=1
else
    CPUS=${SLURM_CPUS_ON_NODE}
fi

module load parallel

#parallel -j 2 fastqsplitter --max-size 10G -S --prefix bc_R{}_ ../barcoded_R{}.fq.gz ::: $(seq 2)

parallel -j 2 zcat ../barcoded_R{}.fq.gz \| split -d --additional-suffix=.fq -l 40000000 - bc_R{}_ ::: $(seq 2)
parallel -j 4 pigz {} ::: $(ls *.fq)
