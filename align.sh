#!/bin/bash

fastqc_bin=$1
data_dir=$2
out_dir=$3
job_script=scripts/jobs/fastqc.sbatch

export num_files=$(find $data_dir -name '*.fastq.gz' | tr ' ' '\n' | wc -l)
# round correctly for truncated arithmetic
export num_files=$((( $num_files + 1 ) / 2 ))

printf \
"#!/bin/bash

#SBATCH -A syb111
#SBATCH -N $num_files
#SBATCH -t 6:00:00
#SBATCH -J fastqc
#SBATCH -o ./scripts/alignment/logs/fastqc.%%J.out
#SBATCH -e ./scripts/alignment/logs/fastqc.%%J.err

module load python

echo fastqc_bin $fastqc_bin
echo data_dir $data_dir
echo out_dir $out_dir

srun -n 2 python ./scripts/alignment/mpi_fastqc.py $fastqc_bin $data_dir $out_dir" \
> $job_script

sbatch $job_script
