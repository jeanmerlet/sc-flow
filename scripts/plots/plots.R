# Libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(dplyr))


# Function arguments
# metadata: metadata from ./data/metadata the user wants to plot (data.frame)
# xvar: the x variable (factor)
# yvar: the y variable to create the violins with (numeric)
# xlab: the label for the x-axis (character)
# ylab: the label for the y-axis (character)
# split_by: create a plot for each level of a factor variable (character)
# condition: whether the user wants to color the violin plots with condition variable (logical)
# mito_cutoff: the mitochondrial cutoff used for the preprocessing (numeric)
# plot_dir: the plotting directory (character)
# plot_name: the name of the plot (character)
# width: width of saved plot (numeric)
# height: height of saved plot (numeric)


# Violin plot

violin_plot <- function(metadata,xvar,yvar,xlab,ylab,color_by,mito_cutoff,plot_dir,plot_name,width,height) {

	metadata[[xvar]] <- as.factor(metadata[[xvar]])
	metadata[[yvar]] <- as.numeric(metadata[[yvar]])
	if (color_by != 'none') {
		plot <- ggplot(data = metadata,aes(x = !!sym(xvar),y = !!sym(yvar),color = !!sym(color_by),fill = !!sym(color_by))) + scale_fill_discrete(name = color_by) + scale_color_discrete(name = color_by)
	} else {
		plot <- ggplot(data = metadata, aes(x = !!sym(xvar), y = !!sym(yvar)))
	}
	plot <- plot +
		geom_violin() +
		geom_hline(yintercept = mito_cutoff, color = "red") +
		coord_flip() +
		xlab(ylab) + 
		ylab(xlab) +
		theme(legend.position = "bottom")
		#theme(strip.text.y = element_text(angle = 0))
	ggsave(paste0(plot_dir,plot_name),plot = plot,width = width,height = height)
	
}


# Bar plot

bar_plot <- function(metadata,xvar,xlab,ylab,color_by,split_by,plot_dir,plot_name,width,height) {
	metadata[[xvar]] <- as.factor(metadata[[xvar]])
	if(color_by != 'none') {
		plot <- ggplot(data = metadata,aes(x = !!sym(xvar),fill = !!sym(color_by))) + scale_fill_discrete(name = color_by)
	} else {
		plot <- ggplot(data = metadata,aes(x = !!sym(xvar)))
	}
	if(split_by != 'none') {
		plot <- plot + facet_grid(as.formula(paste("~",split_by)),scales = "free")
	}
	plot <- plot +
		geom_bar() +
		scale_y_continuous(labels = comma) +
		geom_text(stat = 'count',aes(label = ..count..),vjust = 1.5,size = 4,color = "white") +
		theme(
			strip.text.y = element_text(angle = 0),
			axis.text.y = element_blank(),
			axis.ticks.y = element_blank()
		) +
		xlab(xlab) +
        	ylab(ylab)
	ggsave(paste0(plot_dir,plot_name),plot = plot,width = width,height = height)

}


# Density plot

density_plot <- function(metadata,xvar,xlab,ylab,color_by,split_by,plot_dir,plot_name,width,height) {
	if(color_by != 'none') {
		plot <- ggplot(data = metadata,aes(x = !!sym(xvar),color = !!sym(color_by))) + scale_color_discrete(name = color_by)
	} else {
		plot <- ggplot(data = metadata,aes(x = !!sym(xvar)))
	}
	plot <- plot + geom_density() + scale_x_log10(labels = comma)
	if(split_by != 'none') {
		plot <- plot + facet_wrap(as.formula(paste("~",split_by)))
	}
	plot <- plot +
		theme(strip.text.y = element_text(angle = 0)) +
		xlab(xlab) +
		ylab(ylab)
	ggsave(paste0(plot_dir,plot_name),plot = plot,width = width,height = height)
}


# Qc plot workflow

plot_qc <- function(meta_dir,mito_cutoff,plot_dir,width,height) {
    xvar <- 'sample_ids'
    yvar <- 'Mito'
    meta_path <- paste0(meta_dir, 'cell_meta_filtered.tsv')
    metadata <- read.table(meta_path, sep='\t')
    plot_dir <- paste0(plot_dir, plot_type, '/')
    xlab <- 'Mitochondrial expression (%)'
    ylab <- 'Sample IDs'
    plot_name <- paste0(plot_type, '_cell-meta-filtered', '_xvar-', yvar, '_yvar-', xvar, '.png')
    print("running violin plot")
    violin_plot(metadata, xvar, yvar, xlab, ylab, color_by, mito_cutoff,plot_dir, plot_name, width, height)
    print("violin plot complete")
    xlab <- 'Sample IDs'
    ylab <- ''
    plot_name <- paste0(plot_type, '_cell-meta-filtered', '_xvar-', xvar, '_yvar-total-cells.png')
    print("running bar plot")
    bar_plot(metadata, xvar, xlab, ylab, color_by, split_by,plot_dir, plot_name, width, height)
    print("bar chart complete")
    split_by <- 'sample_ids'
    xvar <- 'nCount_RNA'
    xlab <- 'Total UMIs per cell'
    ylab <- 'Density'
    plot_name <- paste0(plot_type, '_cell-meta-filtered', '_xvar-', xvar, '_yvar-density.png')
    print("running density plot (total UMI)")
    density_plot(metadata,xvar,xlab,ylab,color_by,split_by,plot_dir,plot_name,width,height)
    print("density plot (total UMI) complete")
    xvar <- 'nFeature_RNA'
    xlab <- 'Unique Genes per cell'
    ylab <- 'Density'
    print("running density plot (total genes)")
    plot_name <- paste0(plot_type, '_cell-meta-filtered', '_xvar-', xvar, '_yvar-density.png')
    density_plot(metadata,xvar,xlab,ylab,color_by,split_by,plot_dir,plot_name,width,height)
    print("density plot (total genes) complete")
}


# Umap plot workflow and plotting

plot_umap <- function(meta_dir, plot_dir, width, height, pcs, res, integrated, min_dist, nn, color_by, split_by) {
    if (integrated) {
        obj_type <- 'integrated'
    } else {
        obj_type <- 'preprocessed'
    }
    umap_coords_path <- paste0(meta_dir, 'umap_pcs-', pcs, '_min-dist-', min_dist, '_nn-', nn, '.tsv')
    umap_coords <- read.table(umap_coords_path, sep='\t', header=TRUE)
    if (color_by == 'clusters' | split_by == 'clusters') {
        clusters_path <- paste0(meta_dir, 'clusters_obj-type-', obj_type, '_resolution-', res, '_num-pcs-', pcs, '.tsv')
        clusters <- read.table(clusters_path, sep='\t', header=TRUE)
        umap_coords$clusters <- factor(clusters$clusters, levels=sort(unique(clusters$clusters)))
    }

    # always load sample ids and condition metadata
    metadata_path <- paste0(meta_dir, 'cell_meta_filtered.tsv')
    metadata <- read.table(metadata_path, sep='\t', header=TRUE)
    umap_coords$sample_ids <- metadata$sample_ids
    umap_coords$condition <- metadata$condition
    if (color_by == 'cluster_labels') {
        umap_coords$cluster_labels <- metadata$cluster_labels
    }
    if (color_by == 'none') {
        plot <- ggplot(data = umap_coords, aes(x = UMAP_1, y = UMAP_2))
        color_name <- ''
    } else {
        plot <- ggplot(data = umap_coords, aes(x = UMAP_1, y = UMAP_2, color = !!sym(color_by)))
        color_name <- paste0('_color-by-', color_by)
    }
    plot <- plot + geom_point(alpha=1, size=0.6)+
        guides(color = guide_legend(override.aes = list(title = "",alpha = 1,size = 4)))+
        theme(
            strip.text.y = element_text(angle = 0),
            legend.title = element_blank()
        )+
        xlab("UMAP 1")+
        ylab("UMAP 2")
    if (split_by != 'none') {
        plot <- plot + facet_wrap(as.formula(paste("~", split_by)))
        split_name <- paste0('_split-by-', split_by, '_')
    } else {
        split_name <- '_'
    }
    plot_name <- paste0(
        'umap',color_name,split_name,
        'obj-type-',obj_type,
        '_pcs-',pcs,
        '_res-', res,
        '_min-dist-',min_dist,
        '_nn-',nn,
        '.png')
	ggsave(paste0(plot_dir, 'cells/', plot_name), plot = plot,width = width,height = height)
}


# DEG volcano plot
plot_volcano <- function(de_dir, plot_dir, diff_type, p_value, p_cutoff, log2fc, top_n_genes,width, height) {
    degs_path <- paste0(de_dir, 'degs_diff-type-', diff_type, '_p-value-', p_value, '.tsv')
    degs <- read.table(degs_path, sep='\t', header=TRUE, row.names=1)
    # add significance column
    degs <- degs %>% mutate(
        significance = ifelse(
            avg_log2FC >= log2fc & p_val_adj < p_cutoff,
            "up-regulated",
            ifelse(
                avg_log2FC <= -log2fc & p_val_adj < p_cutoff,
                "down-regulated",
                "non-significant"
            )
        )
    ) %>% data.frame()
    # add labels for most significant genes in each lineage
    labels <- do.call("rbind",lapply(unique(degs[["clusters"]]),function(x){
        up <- degs[degs[["clusters"]] == x & degs$significance == "up-regulated",]
        up <- up[order(up$avg_log2FC,up$p_val_adj,decreasing = c(TRUE,FALSE)),][1:top_n_genes,]
        down <- degs[degs[["clusters"]] == x & degs$significance == "down-regulated",]
        down <- down[order(down$avg_log2FC,down$p_val_adj,decreasing = c(FALSE,FALSE)),][1:top_n_genes,]
        return(rbind(up,down))
    }))
    plot <- ggplot(data = degs,aes(x = avg_log2FC,y = -log10(p_val_adj),col = significance)) +
        geom_point(size = 0.8) +
        facet_wrap(as.formula(paste("~clusters")), scales = "free") +
        scale_color_manual(values = c("down-regulated" = "blue","non-significant" = "black","up-regulated" = "red")) +
        geom_vline(xintercept = c(-log2fc, log2fc),col = "black",linetype = "dashed") +
        geom_hline(yintercept = -log10(p_cutoff),col = "black",linetype = "dashed") +
        geom_label_repel(data = labels,aes(x = avg_log2FC,y = -log10(p_val_adj),label = gene,col = significance),
                         min.segment.length = 0,size = 2,show.legend = FALSE, max.overlaps=100)+
        theme(legend.title = element_blank())+
        guides(color = guide_legend(override.aes = list(alpha = 1,size = 4)))+
        labs(x = expression(log[2](FC))) +
        labs(y = expression(-log[10](adj. ~ p-value)))
    plot_name <- paste0(
        'volcano_diff-type-',diff_type,
        '_min-l2fc-',log2fc,
        '_p-value-',p_value,
        '_p-cutoff-',p_cutoff,
        '_top-',top_n_genes,
        '-genes.png')
	ggsave(paste0(plot_dir, 'de/', plot_name), plot = plot,width = width,height = height)
    warnings()
}
