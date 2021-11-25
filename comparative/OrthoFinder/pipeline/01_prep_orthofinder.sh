#!/usr/bin/bash -l
#SBATCH -p short -C xeon -N 1 -n 96 --mem 24gb --out logs/orthofinder.%A.log

mkdir -p logs
module load orthofinder/2.5.2
opt="" # could change to "-C xeon" and will run on the xeon nodes; # could change this to empty and will run on any node
JOBS=orthofinder_steps.mmseqs.sh
LOG=orthofinder_steps.log
CHUNK=50
export TMPDIR=/scratch
if [ ! -f $LOG ]; then
	orthofinder -op -t 96 -a 96 -f input -S mmseqs -o OrthoFinder_MMSeq > $LOG
	grep ^mmseqs $LOG | grep -v 'commands that must be run' | perl -p -e 's#tmp/#scratch/#g; s/threads 1/threads 8/g'> $JOBS
fi
t=$(wc -l $JOBS | awk '{print $1}')
MAX=$(expr $t / $CHUNK)
echo "t is $t MAX is $MAX"
exit
for n in $(seq $MAX)
do
	START=$(perl -e "printf('%d',1 + $CHUNK * ($n - 1))")
	END=$(perl -e "printf('%d',$CHUNK* $n)")
#	echo "$START,$END for $n"
	run=$(sed -n ${START},${END}p $JOBS)
#		echo "sbatch -p short -N 1 -n 1 --mem 4gb --wrap \"module load orthofinder/2.5.2; $run\""
	sbatch $opt --out logs/MMSeq.$n.log -J MMSeq$n -p short -N 1 -n 8 --mem 4gb --wrap "module load orthofinder/2.5.2; $run"
done
