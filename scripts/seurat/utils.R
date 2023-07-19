suppressPackageStartupMessages({
    library(Seurat)
})

load_seurat_obj <- function(path) {
    obj <- readRDS(file=path)
    return(obj)
}


# add a tsv with bcs and metadata cols to seurat obj metadata
add_metadata <- function(obj, meta_path) {
    metadata <- read.table(meta_path, sep='\t', header=TRUE, row.names=1)
    obj <- AddMetaData(obj, metadata, col.name=NULL)
    return(obj)
}


# Find variable features, scale, PCA
apply_dim_reduce <- function(obj, num_pcs) {
    start_time <- Sys.time()
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj, features=VariableFeatures(object=obj), npcs=num_pcs)
    print(paste0('PCA analysis took ', difftime(Sys.time(), start_time, units='secs')))
    return(obj)
}

