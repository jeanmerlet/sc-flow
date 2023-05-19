#!/bin/bash

# NOTE: run in project root

## set up directory structure
# data
mkdir -p data
cd data
mkdir -p bam
mkdir -p bam/logs
mkdir -p count-matrices
mkdir -p de
mkdir -p metadata
mkdir -p raw
mkdir -p seurat-objects
# plots
cd ..
mkdir -p plots
cd plots
mkdir -p cells
mkdir -p de
mkdir -p genes
mkdir -p qc
mkdir -p qc/reads
# qc
cd ..
mkdir -p qc
# scripts
cd ..
mkdir -p scripts
cd scripts
mkdir -p alignment
mkdir -p download
mkdir -p plots
mkdir -p preprocess
mkdir -p seurat
for dir in ./*/; do mkdir -p -- "$dir/logs"; done
mkdir -p jobs
