#!/bin/bash
#SBATCH --ntasks 48 --nodes 1 --mem 128G -p intel --time 36:00:00 --out iprscan.%A.log

module load iprscan
module load funannotate/git-live
funannotate iprscan -i funannot -m local -c 48 --iprscan_path $(which interproscan.sh)
