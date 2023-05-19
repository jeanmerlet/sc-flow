#!/bin/bash

options=$1

preprocess_job_script=./scripts/jobs/preprocess.sbatch

# preprocess
printf \
"#!/bin/bash

#SBATCH -A syb111
#SBATCH -N 1
#SBATCH -t 6:00:00
#SBATCH -J preprocess
#SBATCH -o ./scripts/preprocess/logs/preprocess.%%J.out
#SBATCH -e ./scripts/preprocess/logs/preprocess.%%J.err

source activate /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/deseq2_andes

srun -n 1 Rscript ./scripts/preprocess/run_preprocess.R $options" \
> $preprocess_job_script

sbatch $preprocess_job_script
