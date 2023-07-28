suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
})


run_diff_exp <- function(obj, de_dir, diff_type, p_value) {
    start_time <- Sys.time()
    if (diff_type == 'cluster') {
        Idents(obj) <- obj@meta.data$clusters
    } else if (diff_type == 'placeholder') {
    }
    obj <- FindVariableFeatures(obj)
    markers <- FindAllMarkers(obj, assay='RNA')
    markers <- markers %>% filter(p_val_adj < p_value) %>% arrange(cluster, desc(avg_log2FC)) %>% data.frame()
    filename <- paste0('degs_diff-type-', diff_type, '_p-value-', p_value, '.tsv')
    out_path <- paste0(de_dir, filename)
    write.table(markers, out_path, sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
    print(paste0('Differential expression runtime: ', difftime(Sys.time(), start_time, unit='mins')))
}
