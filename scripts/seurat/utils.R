library(Seurat)

load_seurat_obj <- function(path) {
    obj <- readRDS(file=path)
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

