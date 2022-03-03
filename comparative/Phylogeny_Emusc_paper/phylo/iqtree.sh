#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 14 --mem 96gb -p short -C xeon  --out fasttree_run.%A.log

#--time 24:00:00 -p intel --out fasttree_run.%A.log

CPU=14
module load IQ-TREE/2.1.3
module unload perl
module unload miniconda2
module load miniconda3
NUM=$(wc -l ../prefix.tab | awk '{print $1}')
source ../config.txt

ALN=../$PREFIX.${NUM}_taxa.$HMM.aa.fasaln
iqtree2 -nt $CPU -s $ALN --prefix $PREFIX.${NUM}_taxa.$HMM
TREE1=$PREFIX.${NUM}_taxa.$HMM.iqtree
TREE2=$PREFIX.${NUM}_taxa.$HMM.long_iq.tre
if [ -s $TREE1 ]; then
	perl ../PHYling_unified/util/rename_tree_nodes.pl $TREE1 ../prefix.tab > $TREE2
fi
