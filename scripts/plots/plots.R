# Libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(scales))


# Function arguments
# metadata: metadata from ./data/metadata the user wants to plot (data.frame)
# xvar: the x variable (factor)
# yvar: the y variable to create the violins with (numeric)
# xlab: the label for the x-axis (character)
# ylab: the label for the y-axis (character)
# facet_var: create a plot for each level of a factor variable (character)
# condition: whether the user wants to color the violin plots with condition variable (logical)
# mito_cutoff: the mitochondrial cutoff used for the preprocessing (numeric)
# plot_dir: the plotting directory (character)
# plot_name: the name of the plot (character)
# width: width of saved plot (numeric)
# height: height of saved plot (numeric)


# Violin plot

violin_plot <- function(metadata,xvar,yvar,xlab,ylab,condition,mito_cutoff,plot_dir,plot_name,width,height) {

	metadata[[xvar]] <- as.factor(metadata[[xvar]])
	metadata[[yvar]] <- as.numeric(metadata[[yvar]])
	if(condition) {
		plot <- plot + ggplot(data = metadata,aes_string(x = xvar,y = yvar,color = condition,fill = condition)) + scale_fill_discrete(name = "Condition") + scale_color_discrete(name = "Condition")
	} else {
        
        plot <- ggplot(data=metadata, aes(x=!!sym(xvar), y=!!sym(yvar)))
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

bar_plot <- function(metadata,xvar,xlab,ylab,condition,plot_dir,plot_name,width,height) {
	metadata[[xvar]] <- as.factor(metadata[[xvar]])
	if(condition) {
		plot <- plot + ggplot(data = metadata,aes(x = xvar,fill = condition)) + facet_grid(.~condition,scales = "free") + scale_fill_discrete(name = "Condition")
	} else {
        plot <- ggplot(data=metadata, aes(x=!!sym(xvar)))
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

density_plot <- function(metadata,xvar,xlab,ylab,facet_var,condition,plot_dir,plot_name,width,height) {
	if(condition) {
		plot <- plot + ggplot(data = metadata,aes(x = xvar,color = condition)) + scale_color_discrete(name = "Condition")
	} else {
        plot <- ggplot(data=metadata, aes(x=!!sym(xvar)))
	}
	plot <- plot +
		geom_density() +
		scale_x_log10(labels = comma) +
        facet_wrap(as.formula(paste("~", facet_var))) +
		#facet_wrap(!!sym(facet_var)~.) +
		theme(strip.text.y = element_text(angle = 0)) +
		xlab(xlab) +
		ylab(ylab)
	ggsave(paste0(plot_dir,plot_name),plot = plot,width = width,height = height)
}


# Qc plot workflow

plot_qc <- function(metadata, xvar, yvar, condition, mito_cutoff,
                    plot_dir, plot_name, width, height) {
    xvar <- 'sample_ids'
    yvar <- 'Mito'
    meta_path <- paste0(meta_dir, 'cell_meta_filtered.tsv')
    metadata <- read.table(meta_path, sep='\t')
    plot_dir <- paste0(plot_dir, plot_type, '/')
    xlab <- 'Mitochondrial expression (%)'
    ylab <- 'Sample IDs'
    plot_name <- paste0(plot_type, '_cell-meta-filtered', '_xvar-', yvar, '_yvar-', xvar, '.png')
    violin_plot(metadata, xvar, yvar, xlab, ylab, condition, mito_cutoff,
                plot_dir, plot_name, width, height)
    xlab <- 'Sample IDs'
    ylab <- ''
    plot_name <- paste0(plot_type, '_cell-meta-filtered', '_xvar-', xvar, '_yvar-total-cells.png')
    bar_plot(metadata, xvar, xlab, ylab, condition, plot_dir, plot_name, width, height)
    facet_var <- 'sample_ids'
    xvar <- 'nCount_RNA'
    xlab <- 'Total UMIs per cell'
    ylab <- 'Density'
    plot_name <- paste0(plot_type, '_cell-meta-filtered', '_xvar-', xvar, '_yvar-density.png')
    density_plot(metadata,xvar,xlab,ylab,facet_var,condition,plot_dir,plot_name,width,height)
    xvar <- 'nFeature_RNA'
    xlab <- 'Unique Genes per cell'
    ylab <- 'Density'
    plot_name <- paste0(plot_type, '_cell-meta-filtered', '_xvar-', xvar, '_yvar-density.png')
    density_plot(metadata,xvar,xlab,ylab,facet_var,condition,plot_dir,plot_name,width,height)
}


# Umap plot workflow and plotting

plot_umap <- function(meta_dir, plot_dir, width, height, pcs, res, integrated, min_dist, nn) {
    if (integrated) {
        obj_type <- 'integrated'
    } else {
        obj_type <- 'preprocessed'
    }
    umap_coords_path <- paste0(meta_dir, 'umap_pcs-', pcs, '_min-dist-', min_dist, '_nn-', nn, '.tsv')
    clusters_path <- paste0(meta_dir, 'clusters_obj-type-', obj_type, '_resolution-', res, '_num-pcs-', pcs, '.tsv')
    metadata <- read.table(umap_coords_path, sep='\t', header=TRUE)
    clusters <- read.table(clusters_path, sep='\t', header=TRUE)
    metadata$clusters <- factor(clusters$clusters, levels=sort(unique(clusters$clusters)))
    plot <- ggplot(data = metadata, aes(x = UMAP_1, y = UMAP_2, color = clusters))+
        geom_point(alpha = 1, size = 0.8)+
        guides(color = guide_legend(override.aes = list(title = "",alpha = 1,size = 4)))+
        theme(
            strip.text.y = element_text(angle = 0),
            legend.title = element_blank()
        )+
        xlab("UMAP 1")+
        ylab("UMAP 2")
    plot_name <- paste0('umap-clusters_obj-type-', obj_type, '_pcs-', pcs, '_res-', res,
                        '_min-dist-', min_dist, '_nn-', nn, '.png')
	ggsave(paste0(plot_dir, 'cells/', plot_name), plot = plot,width = width,height = height)
}
