


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


add_cluster_labels_new <- function(labels_path, cell_meta_path, clusters_path) {
    # read labels
    labels <- read.table(labels_path, sep='\t', header=TRUE, row.names=1)
    # read cell_meta
    cell_meta <- read.table(cell_meta_path, sep='\t', header=TRUE, row.names=1)
    # read clusters
    clusters <- read.table(clusters_path, sep='\t', header=TRUE, row.names=1)
    # add labels to clusters
    for (bc in rownames(clusters)) {
        cluster <- as.character(clusters[bc, 'clusters'])
        cluster_label <- labels[cluster,]
        clusters[bc, 'labeled_clusters'] <- cluster_label
    }
    # drop the cluster numbers
    clusters$clusters <- NULL
    colnames(clusters) <- c('clusters')
    # write to file
    out_path <- gsub('clusters_obj-type', 'annotated_clusters_obj-type', clusters_path)
    write.table(clusters, out_path, sep='\t', quote=FALSE)
}
