import subprocess
import sys, os, re

# script arguments
star_bin = sys.argv[1]
out_dir = sys.argv[2]
genome_fasta = sys.argv[3]
genome_annot = sys.argv[4]

# alignment function
def star_index(star_bin, out_dir, genome_fasta, genome_annot):
    subprocess.run([star_bin,
                    '--runThreadN', '8',
                    '--runMode', 'genomeGenerate',
                    '--genomeDir', out_dir,
                    '--genomeFastaFiles', genome_fasta,
                    '--sjdbGTFfile', genome_annot])

star_index(star_bin, out_dir, genome_fasta, genome_annot)
