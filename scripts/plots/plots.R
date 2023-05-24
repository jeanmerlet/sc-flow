
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


# Violin plot function

violin_plot <- function(metadata,xvar,yvar,xlab,ylab,condition,mito_cutoff,plot_dir,plot_name,width,height) {

	xvar <- enquo(xvar)
	yvar <- enquo(yvar)
	if(condition) {
		plot <- plot + ggplot(data = metadata,aes_string(x = xvar,y = yvar,color = condition,fill = condition)) + scale_fill_discrete(name = "Condition") + scale_color_discrete(name = "Condition")
	} else {
		plot <- ggplot(data = metadata,aes(x = !!xvar,y = !!yvar))
	}
	plot <- plot +
		geom_violin() +
		geom_hline(yintercept = mito_cutoff,color = "red") +
		coord_flip() +
		xlab(xlab) + 
		ylab(ylab) +
		theme(legend.position = "bottom")
		#theme(strip.text.y = element_text(angle = 0))
	ggsave(paste0(plot_dir,plot_name),plot = plot,width = width,height = height)
	
}


# Bar plot function

bar_plot <- function(metadata,xvar,xlab,condition,plot_dir,plot_name,width,height) {

	plot <- ggplot()
	if(condition) {
		plot <- plot + ggplot(data = metadata,aes(x = xvar,fill = condition)) + facet_grid(.~condition,scales = "free") + scale_fill_discrete(name = "Condition")
	} else {
		plot <- plot + ggplot(data = metadata,aes(x = xvar,fill = condition))
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


# Density plot function

density_plot <- function(metadata,xvar,xlab,ylab,facet_var,condition,plot_dir,plot_name,width,height) {

	plot <- ggplot()
	if(condition) {
		plot <- plot + ggplot(data = metadata,aes(x = xvar,color = condition)) + scale_color_discrete(name = "Condition")
	} else {
		plot <- plot + ggplot(data = metadata,aes(x = xvar))
	}
	plot <- plot +
		geom_density() +
		scale_x_log10(labels = comma) +
		facet_wrap(facet_var~.) +
		theme(strip.text.y = element_text(angle = 0)) +
		xlab(xlab) +
		ylab(ylab)
	ggsave(paste0(plot_dir,plot_name),plot = plot,width = width,height = height)

}
