
# Libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(scales))


violin_plot <- function(metadata,xvar,yvar,xlab,ylab,condition,mito_cutoff,plot_dir,plot_name,width,height) {

	plot <- ggplot()
	if(condition) {
		plot <- plot + ggplot(data = metadata,aes_string(x = xvar,y = yvar,color = condition,fill = condition)) +
			theme(legend.position = "bottom")
	} else {
		plot <- plot + ggplot(data = metadata,aes(x = xvar,y = yvar))
	}
	plot <- plot +
		geom_violin() +
		geom_hline(yintercept = mito_cutoff,color = "red") +
		coord_flip()+
		xlab(xlab)+
		ylab(ylab)
		#theme(strip.text.y = element_text(angle = 0))
	ggsave(paste0(plot_dir,plot_name),plot = plot,width = width,height = height)
	
}


#apply_qc_plots_v2 <- function(type,condition) {
#
#	# type: "raw" or "preprocessed"
#	# condition: TRUE or FALSE (default: FALSE)
#	obj <- readRDS(",/data/seurat-objects/",type,".rds")
#	metadata <- obj@meta.data
#	metadata$sample_id <- factor(metadata$sample_id,levels = sort(unique(metadata$sample_id)),labels = sort(unique(metadata$sample_id)))
#	if(condition) {
#		metadata$condition <- factor(metadata$condition,levels = sort(unique(metadata$condition)),labels = sort(unique(metadata$condition)))
#	}
#	# violin plots for mitochondrial content
#	# number of cells bar chart
#	plot <- ggplot()
#	if(condition) {
#		plot <- plot + ggplot(data = metadata,aes(x = sample_id,fill = condition)) + facet_grid(.~condition,scales = "free")
#	} else {
#		plot <- plot + ggplot(data = metadata,aes(x = sample_id,fill = condition))
#	}
#	plot <- plot +
#		geom_bar() +
#		scale_y_continuous(labels = comma) +
#		geom_text(stat = 'count',aes(label = ..count..),vjust = 1.5,size = 4,color = "white") +
#		xlab("") +
#		scale_fill_discrete(name = "") +
#		theme(
#			strip.text.y = element_text(angle = 0),
#			axis.text.y = element_blank(),
#			axis.ticks.y = element_blank()
#		) +
#		ylab("Total cells")
#	ggsave(paste0("./plots/qc/",type,"_total_cell_numbers.png"),plot = plot,width = 10,height = 5)
#	# distribution of total umi per cell
#	if(type == "raw") {
#		metadata <- metadata[metadata$nCount_RNA > 10,]
#	}
#	plot <- ggplot()
#	if(condition) {
#		plot <- plot + ggplot(data = metadata,aes(x = nCount_RNA,color = condition))
#	} else {
#		plot <- plot + ggplot(data = metadata,aes(x = nCount_RNA))
#	}
#	plot <- plot +
#		geom_density() +
#		scale_x_log10(labels = comma) +
#		facet_wrap(sample_id~.) +
#		theme(strip.text.y = element_text(angle = 0)) +
#		xlab("UMI count") +
#		ylab("Density")
#	ggsave(paste0("./plots/qc/",type,"_total_umi_distribution.png"),plot = plot,width = 10,height = 6)
#	# distribution of total genes per cell
#	plot <- ggplot()
#	if(condition) {
#		plot <- plot + ggplot(data = metadata,aes(x = nFeature_RNA,color = condition))
#	} else {
#		plot <- plot + ggplot(data = metadata,aes(x = nFeature_RNA))
#	}
#	plot <- plot +
#		geom_density() +
#		scale_x_log10(labels = comma) +
#		facet_wrap(sample_id~.) +
#		theme(strip.text.y = element_text(angle = 0)) +
#		xlab("Gene count") +
#		ylab("Density")
#	ggsave(paste0("./plots/qc/",type,"_total_gene_distribution.png"),plot = plot,width = 10,height = 6)
#
#}
#
#apply_qc_plots_v1(type = opt$type,condition = opt$condition)

