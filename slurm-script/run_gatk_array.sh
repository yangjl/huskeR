#!/bin/bash -l
#SBATCH -D /Users/yangjl/Documents/Github/huskeR
#SBATCH -o /Users/yangjl/Documents/Github/huskeR/slurm-log/testout-%j.txt
#SBATCH -e /Users/yangjl/Documents/Github/huskeR/slurm-log/err-%j.txt
#SBATCH -J runarray
#SBATCH --array=1-1
#SBATCH --mail-user=
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL #email if fails
set -e
set -u

module load java
module load bwa
sh slurm-script/run_gatk_1.sh
