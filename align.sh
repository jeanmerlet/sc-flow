#!/bin/bash

fastqc_bin=$1
star_bin=$2
genome_dir=$3
bc_whitelist=$4
umi_len=$5
data_dir=./data/raw
out_fastqc_dir=./qc/reads
out_align_dir=./data/bam
out_alignqc_dir=./qc/alignment


export num_files=$(find $data_dir -name '*.fastq.gz' | tr ' ' '\n' | wc -l)
# round correctly for truncated arithmetic
export num_paired_files=$((( $num_files + 1 ) / 2 ))

# fastqc
printf \
"#!/bin/bash

#SBATCH -A syb111
#SBATCH -N $num_paired_files
#SBATCH -t 6:00:00
#SBATCH -J fastqc
#SBATCH -o ./scripts/alignment/logs/fastqc.%%J.out
#SBATCH -e ./scripts/alignment/logs/fastqc.%%J.err

module load python

echo fastqc_bin $fastqc_bin

srun -n 2 python ./scripts/alignment/mpi_fastqc.py $fastqc_bin $data_dir $out_fastqc_dir" \
#rm $out_fastqc_dir/*.zip
#multiqc $out_fastqc_dir
#mv $out_fastqc_dir/multiqc_report.html $out_fastqc_dir/fastqc_report.html" \
> $fastqc_job_script

sbatch $fastqc_job_script


# alignment and alignment qc
printf \
"#!/bin/bash

#SBATCH -A SYB111
#SBATCH -N $num_paired_files
#SBATCH -t 6:00:00
#SBATCH -J STAR_align
#SBATCH -o ./logs/align.%%J.out
#SBATCH -e ./logs/align.%%J.err

module load python

echo star_bin $star_bin

srun -n 1 python ./scripts/alignment/mpi_align.py $star_bin $genome_dir $bc_whitelist $data_dir $out_align_dir $umi_len

source load /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/python_andes
srun -N 1 -n 1 multiqc $out_align_dir
mv $out_align_dir/multiqc_report.html $out_alignqc_dir/align_report.html" \
> $align_job_script

#sbatch $align_job_script
