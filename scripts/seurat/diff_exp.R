suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(future)
    library(data.table)
    library(doParallel)
})

options(future.globals.maxSize = 100000 * 1024^2,future.rng.onMisue = "ignore")

run_diff_exp <- function(obj, res, de_dir, diff_type, condition_group, p_value, use_annotated) {
    start_time <- Sys.time()
    Idents(obj) <- obj@meta.data$clusters
    obj <- FindVariableFeatures(obj)
    #TODO: parallelize clusters by cluster
    if (diff_type == 'cluster') {
        plan('multicore',workers = 8)
        markers <- FindAllMarkers(obj, assay='RNA')
        names(markers)[which(names(markers) == "cluster")] <- "clusters"
    } else if (diff_type == 'condition') {
        cl <- makeCluster(8)
        registerDoParallel(cl)
        markers <- foreach(clust = sort(unique(obj$clusters)),.packages = c("Seurat")) %dopar% {
            cond1 <- unique(obj@meta.data$condition)[1]
            cond2 <- unique(obj@meta.data$condition)[2]
            if (sum(obj[,obj@meta.data$clusters == clust]@meta.data$condition == cond1) > 0 &
              sum(obj[,obj@meta.data$clusters == clust]@meta.data$condition == cond2) > 0) {
                cluster_markers <- FindMarkers(obj, ident.1=condition_group, group.by="condition", subset.ident=clust)
                cluster_markers$gene <- rownames(cluster_markers)
                cluster_markers$clusters <- clust
                return(cluster_markers)
            } else {
                print(paste0('cluster ', clust, ' has no cells in one of the conditions, excluding it from results'))
                console.flush()
                return(NULL)
            }
        }
        stopCluster(cl)
        markers <- markers[!sapply(result, is.null)]
        markers <- rbindlist(markers)
    }
    markers <- markers %>% filter(p_val_adj < p_value) %>% arrange(clusters, desc(avg_log2FC)) %>% data.frame()
    if (use_annotated) {
        annot_name <- 'annotated_'
    } else {
        annot_name <- ''
    }
    filename <- paste0(annot_name, 'degs_diff-type-', diff_type, '_p-value-', p_value, '_res-', res, '.tsv')
    out_path <- paste0(de_dir, filename)
    # row.names false because there is a gene col
    write.table(markers, out_path, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
    print(paste0('Differential expression runtime: ', difftime(Sys.time(), start_time, unit='mins')))
}

diff_exp <- function(meta_dir,resolution, pcs, use_integrated, de_dir,
                     diff_type, condition_group, p_value, use_annotated) {
    # load imputed obj for generating DEGs
    obj_type <- 'imputed'
    obj_path <- paste0(obj_dir,obj_type,'.rds')
    obj <- load_seurat_obj(obj_path)
    # load cluster metadata from either preprocessed or integrated obj
    if (use_integrated) {
        obj_type <- 'integrated'
    } else {
        obj_type <- 'preprocessed'
    }
    # use cluster numbers or cluster labels
    if (use_annotated) {
        meta_path <- paste0(meta_dir, 'annotated_clusters_obj-type-', obj_type, '_resolution-',
                            resolution,'_num-pcs-', pcs, '.tsv')
    } else {
        meta_path <- paste0(meta_dir, 'clusters_obj-type-', obj_type, '_resolution-',
                            resolution,'_num-pcs-', pcs, '.tsv')
    }
    obj <- add_metadata(obj, meta_path)
    run_diff_exp(obj, resolution, de_dir, diff_type, condition_group, p_value, use_annotated)
}
