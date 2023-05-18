#!/bin/bash

fastqc_bin=$1
star_bin=$2
genome_dir=$3
bc_whitelist=$4
umi_len=$5

data_dir=./data/raw
out_align_dir=./data/bam
out_qc_dir=./qc
out_mtx_dir=./data/count-matrices
fastqc_job_script=./scripts/jobs/fastqc.sbatch
align_job_script=./scripts/jobs/align.sbatch


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
echo fastqc_bin: $fastqc_bin
echo
srun -n $num_files python ./scripts/alignment/mpi_fastqc.py $fastqc_bin $data_dir $out_fastqc_dir

source activate /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/python_andes
rm $out_fastqc_dir/*.html
srun -N 1 -n 1 multiqc $out_qc_dir -n fastqc_report.html -o $out_qc_dir --no-data-dir
rm $out_fastqc_dir/*.zip" \
> $fastqc_job_script

sbatch $fastqc_job_script


# alignment and alignment qc
printf \
"#!/bin/bash

#SBATCH -A SYB111
#SBATCH -N $num_paired_files
#SBATCH -t 6:00:00
#SBATCH -J STAR_align
#SBATCH -o ./scripts/alignment/logs/align.%%J.out
#SBATCH -e ./scripts/alignment/logs/align.%%J.err

module load python
echo star_bin: $star_bin
echo
srun -n $num_paired_files python ./scripts/alignment/mpi_align.py $star_bin $genome_dir $bc_whitelist $data_dir $out_align_dir $umi_len

source activate /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/python_andes
srun -N 1 -n 1 multiqc $out_qc_dir -n align_report.html -o $out_qc_dir --no-data-dir
mv $out_align_dir/*Log*.out $out_align_dir/logs
mv $out_align_dir/*Solo.out $out_mtx_dir" \
> $align_job_script

sbatch $align_job_script
