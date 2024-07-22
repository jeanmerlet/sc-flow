# sc-flow

## setup
This package is designed for use on frontier and assumes you have a working conda or miniconda installation.
To generate your own copy of our R env, do the following steps:
1. Run ```git clone https://github.com/jeanmerlet/sc-flow.git```
2. In /env, run ```echo "prefix: /path/to/env/dir" >> frontier_seurat.yml```
3. In /env, run ```conda env create -f frontier_seurat.yml```
4. In the root directory, run ```./project_init.sh```
5. Place your fastqs in /data/raw
6. Create and fill a /data/metadata/sample_meta.tsv file (see [metadata](#metadata))

## Broad order of steps
In general, a typical scRNA-seq workflow will go through these steps (in order).
1. align
2. preprocess
3. cluster
4. plot
5. diff\_exp

Each of these corresponds to an argument for the ```--workflow``` command in the sc-flow pipeline.
1. align: ```--workflow align```
2. preprocess: ```--workflow preprocess```
3. cluster: ```--workflow cluster```
4. plot: ```--workflow plot```
5. diff\_exp: ```--workflow diff_exp```

## Example workflow
1. ```Rscript sc-flow.R --workflow align``` (options that we have forgotten)
2. ```Rscript sc-flow.R --workflow preprocess --species mouse --load_meta```
3. ```Rscript sc-flow.R --workflow plot --plot_type qc --species mouse --color_by condition```
4. ```Rscript sc-flow.R --workflow plot --plot_type umap --color_by sample_ids```
5. ```Rscript sc-flow.R --workflow plot --plot_type umap --color_by condition```
6. ```Rscript sc-flow.R --workflow cluster --resolution 0.5```
7. ```Rscript sc-flow.R --workflow plot --plot_type umap --color_by clusters```
8. ```Rscript sc-flow.R --workflow plot --plot_type umap --color_by clusters --split_by condition```
9. ```Rscript sc-flow.R --workflow diff_exp --diff_type cluster```
10. ```Rscript sc-flow.R --workflow diff_exp --diff_type condition``` (NOT YET IMPLEMENTED)


## NOTE: talk about interactivity with ```--run\_r```
## ANOTHER NOTE: make sure ```--run\_r``` works with some weird species options stuff for preprocess / qc plot
## NOTE: tell users where there data is located (especially for plots, count matrices, etc.)
