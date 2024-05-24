from mpi4py import MPI
import subprocess
import sys, os, re

# script arguments
star_bin = sys.argv[1]
data_dir = sys.argv[2]
out_dir = sys.argv[3]
genome_dir = sys.argv[4]
bc_whitelist = sys.argv[5]
bc_length = sys.argv[6]
umi_start = sys.argv[7]
umi_len = sys.argv[8]

# alignment function
def sc_star_align(fastq1, fastq2, prefix, genome_dir, bc_whitelist, bc_length, umi_start, umi_len, strandedness):
    out_prefix = os.path.join(out_dir, prefix + '_')
    subprocess.run([star_bin,
                    '--runThreadN', '16',
                    '--genomeDir', genome_dir,
                    '--soloCBwhitelist', bc_whitelist,
                    '--soloCBmatchWLtype', '1MM_multi_Nbase_pseudocounts',
                    '--soloCBlen', bc_length,
                    '--soloUMIstart', umi_start,
                    '--soloUMIlen', umi_len,
                    '--soloBarcodeReadLength', '0',
                    '--soloCellFilter','EmptyDrops_CR',
                    '--soloType CB_UMI_Simple',
                    '--soloFeatures', 'Gene',
                    '--soloUMIdedup', '1MM_CR',
                    '--soloUMIfiltering', 'MultiGeneUMI_CR',
                    '--clipAdapterType', 'CellRanger4',
                    '--outFilterScoreMin', '30',
                    '--readStrand', strandedness,
                    '--readFilesCommand', 'zcat',
                    '--readFilesIn', fastq2, fastq1,
                    '--outFileNamePrefix', out_prefix,
                    '--outSAMtype', 'BAM', 'SortedByCoordinate'])


fastq_paths = os.listdir(data_dir)
fastq_paths.sort()

paired_fastq_list = []
pair = []
for i, fastq_path in enumerate(fastq_paths):
    pair.append(os.path.join(data_dir, fastq_path))
    if i % 2 == 1:
        paired_fastq_list.append(pair)
        pair = []

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, pair in enumerate(paired_fastq_list):
    if i % size != rank: continue
    fastq_path1, fastq_path2 = pair
    head, fastq = os.path.split(fastq_path1)
    prefix = re.search('(.*)[rR][12].*', fastq).groups()[0]
    if prefix[-1] == '_': prefix = prefix[:-1]
    print(f'{prefix} (task number {i}) being aligned by {rank} of {size}')
    sc_star_align(fastq_path1, fastq_path2, prefix, genome_dir, bc_whitelist, bc_length, umi_start, umi_len)
