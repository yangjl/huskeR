#!/bin/bash -l
#SBATCH -D /Users/yangjl/Documents/Github/huskeR/test
#SBATCH -o /Users/yangjl/Documents/Github/huskeR/test/slurm-log/testout-%j.txt
#SBATCH -e /Users/yangjl/Documents/Github/huskeR/test/slurm-log/err-%j.txt
#SBATCH -J runarray
#SBATCH --array=1-1
#SBATCH --mail-user=
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL #email if fails
set -e
set -u

module load java
module load bwa
module load samtools
sh slurm-script/run_test_$SLURM_ARRAY_TASK_ID.sh
