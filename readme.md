# sc-flow

## setup
This package is designed for use on frontier and assumes you have a working conda or miniconda installation.
To generate your own copy of our R env, run the following commands:
1. ```git clone https://github.com/jeanmerlet/sc-flow.git```
2. In /env ```echo "prefix: /path/to/your/env" >> frontier_seurat.yml```
3. ```conda env create -f frontier_seurat.yml```
4. ```./project_init.sh```
5. Place your fastqs in the /data/raw directory
6. Create and fill a /data/metadata/sample_meta.tsv file (see [metadata])

# NOTE: test yml env creation ^^

## Broad order of steps
In general, a typical scRNA-seq workflow will go through these steps (in order).
1. align

## NOTE: fix align arguments / options to appear in argparse

1. preprocess
2. cluster
3. plot
4. diff\_exp

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

## NOTE: fix qc plots to accept color\_by instead of condition and look at split\_by for density plots

1. ```Rscript sc-flow.R --workflow plot --plot_type umap --color_by sample_ids```
2. ```Rscript sc-flow.R --workflow plot --plot_type umap --color_by condition```
3. ```Rscript sc-flow.R --workflow cluster --resolution 0.5```
4. ```Rscript sc-flow.R --workflow plot --plot_type umap --color_by clusters```
5. ```Rscript sc-flow.R --workflow plot --plot_type umap --color_by clusters --split_by condition```
6. ```Rscript sc-flow.R --workflow diff_exp --diff_type cluster```
7. ```Rscript sc-flow.R --workflow diff_exp --diff_type condition``` (NOT YET IMPLEMENTED)

## NOTE: implement diff\_exp by condition, change diff\_type cluster to clusters


## align

## preprocess
The preprocess workflow has the following options.
1. ```rare_gene_cutoff``` with possible arguments any integer between 0 to 100 (corresponding to the %).
2. ```mito_cutoff```
3. ```upper_umi_cutoff```
4. ```load_meta```


## Metadata

rare_gene_cutoff, mito_cutoff, upper_umi_cutoff, load_meta)


## NOTE: talk about interactivity with ```--run\_r```
## ANOTHER NOTE: make sure ```--run\_r``` works with some weird species options stuff for preprocess / qc plot
