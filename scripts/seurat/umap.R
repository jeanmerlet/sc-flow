library(Seurat)


#TODO: move integrate logic to utils
gen_umap <- function(meta_dir, npcs, min_dist, nn, use_integrated) {
    if (use_integrated) {
        obj_path <- './data/seurat-objects/integrated.rds'
    } else {
        obj_path <- './data/seurat-objects/preprocessed.rds'
    }
    obj <- load_seurat_obj(obj_path)
    obj <- apply_dim_reduce(obj, npcs)
    obj <- RunUMAP(obj,reduction = 'pca',verbose = FALSE,dims = 1:npcs, umap.method = 'umap-learn',
                   metric = 'correlation',min.dist = min_dist, n.neighbors = nn)
    umap_data <- as.data.frame(Embeddings(object = obj[['umap']]))
    names(umap_data) <- c('UMAP_1','UMAP_2')
    out_path <- paste0(meta_dir, 'umap_pcs-', npcs, '_min-dist-', min_dist, '_nn-', nn, '.tsv')
    write.table(umap_data, out_path, sep = '\t', row.names = TRUE, col.names = TRUE, quote=FALSE)
}
