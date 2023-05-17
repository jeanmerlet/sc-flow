#!/bin/bash

fastqc_bin=$1
data_dir=$2
out_dir=$3

export num_files=$(find $data_dir -name '*.fastq.gz' | tr ' ' '\n' | wc -l)
# round correctly for truncated arithmetic
export num_files=$( $num_files + 1 / $num_files )

#TODO: log the job script
sbatch \
    -A syb111 \
    -N $num_files \
    -t 6:00:00 \
    -J fastqc \
    -o ./logs/fastqc.%J.out \
    -e ./logs/fastqc.%J.err \
    module load python \
    srun -n 2 python ./mpi_fastqc.py $fastqc_bin $data_dir $out_dir
