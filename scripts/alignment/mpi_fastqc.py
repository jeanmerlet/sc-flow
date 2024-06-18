from mpi4py import MPI
import subprocess
import sys, os

# script arguments
fastqc_bin = sys.argv[1]
data_dir = sys.argv[2]
out_dir = sys.argv[3]

# get list of fastq paths
fastq_paths = [os.path.join(data_dir, path) for path in os.listdir(data_dir)]

# run fastqc command, given a file
def run_fastqc(path):
    subprocess.run([fastqc_bin, '--threads=16', f'--outdir={out_dir}', path])

# distribute over processes
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, fastq_path in enumerate(fastq_paths):
    # distribute across ranks
    if i % size != rank: continue
    head, fastq = os.path.split(fastq_path)
    print(f'running fastqc on {fastq} (task number {i}) by {rank} of {size}')
    # run the command line command on this particular mpi rank
    run_fastqc(fastq_path)
