suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(future)
    library(data.table)
    library(doParallel)
})

options(future.globals.maxSize = 100000 * 1024^2,future.rng.onMisue = "ignore")

run_diff_exp <- function(obj, de_dir, diff_type, condition_group, p_value) {
    start_time <- Sys.time()
    print(start_time)
    Idents(obj) <- obj@meta.data$clusters
    print("changed Idents")
    obj <- FindVariableFeatures(obj)
    print("found variable features")
    if (diff_type == 'cluster') {
        plan('multicore',workers = 8)
        markers <- FindAllMarkers(obj, assay='RNA')
    } else if (diff_type == 'condition') {
	print(diff_type)
        print(condition_group)
	print("starting DE")
        cl <- makeCluster(8)
        registerDoParallel(cl)
        markers <- foreach(clust = sort(unique(obj$clusters)),.packages = c("Seurat")) %dopar% {
            test <- FindMarkers(obj,ident.1 = condition_group,group.by = "condition",subset.ident = clust)
            test$gene <- rownames(test)
            test$clusters <- clust
            return(test)
        }
        stopCluster(cl)
        markers <- rbindlist(markers)
        print("completed DE")	
    }
    markers <- markers %>% filter(p_val_adj < p_value) %>% arrange(clusters, desc(avg_log2FC)) %>% data.frame()
    filename <- paste0('degs_diff-type-', diff_type, '_p-value-', p_value, '.tsv')
    out_path <- paste0(de_dir, filename)
    write.table(markers, out_path, sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
    print(paste0('Differential expression runtime: ', difftime(Sys.time(), start_time, unit='mins')))
}
