


add_cluster_labels <- function(labels_path, cell_meta_path, clusters_path) {
    # read labels
    labels <- read.table(labels_path, sep='\t', header=TRUE, row.names=1)
    # read cell_meta
    cell_meta <- read.table(cell_meta_path, sep='\t', header=TRUE, row.names=1)
    # read clusters
    clusters <- read.table(clusters_path, sep='\t', header=TRUE, row.names=1)
    # add labels to cell_meta
    for (bc in rownames(clusters)) {
        cluster <- as.character(clusters[bc,])
        cluster_label <- labels[cluster,]
        cell_meta[bc, 'cluster_labels'] <- cluster_label
    }
    # write cell_meta
    write.table(cell_meta, cell_meta_path, sep='\t', quote=FALSE)
}
