library(Seurat)

#TODO: subsets and subset names
apply_clustering <- function(obj, out_dir, num_pcs, res, use_integrated) {
    start_time <- Sys.time()
    obj <- FindNeighbors(obj, dims=1:num_pcs)
    obj <- FindClusters(obj, resolution=res)
    if (use_integrated) {
        obj_type <- 'integrated'
    } else {
        obj_type <- 'preprocessed'
    }
    out_path <- paste0(out_dir, 'clusters_obj-type-', obj_type, '_resolution-', res, '_num-pcs-', num_pcs, '.tsv')
    data_frame <- data.frame(obj$seurat_clusters)
    names(data_frame) <- c('clusters')
    write.table(data_frame, row.names=TRUE, col.names=NA, file=out_path, sep='\t', quote=FALSE)
    print(paste0('Clustering took ', difftime(Sys.time(), start_time, units='secs')))
}


# Cluster workflow
cluster <- function(out_dir, pcs, resolution, use_integrated) {
    if (use_integrated) {
        obj_path <- './data/seurat-objects/integrated.rds'
    } else {
        obj_path <- './data/seurat-objects/preprocessed.rds'
    }
    obj <- load_seurat_obj(obj_path)
    obj <- apply_dim_reduce(obj, pcs)
    apply_clustering(obj, out_dir, pcs, resolution, use_integrated)
}

