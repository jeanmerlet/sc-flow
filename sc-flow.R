# Libraries
suppressPackageStartupMessages({
    library(optparse)
    library(future)
})
# /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/deseq2_andes
# /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_seurat

raw_args <- paste0(commandArgs(trailingOnly = TRUE), collapse=' ')

options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = 100 * 1024^3)

mtx_dir <- './data/count-matrices/'
obj_dir <- './data/seurat-objects/'
plot_dir <- './plots/'
meta_dir <- './data/metadata/'


#TODO: standardize timing print statements
# Argument parser
option_list <- list(
    make_option(
        c('--run_r'),
        type='logical',
        action='store_true',
        default=FALSE,
        help='if true, a job script is calling sc-flow to run a workflow'
    ),
    make_option(
        c('--workflow'),
        type='character',
        default=NULL,
        help='what kind of workflow to run'
    ),
    make_option(
        c('--obj_path'),
        type='character',
        default=NULL,
        help='explicit path to a seurat object'
    ),
    make_option(
        c('--out_impute_dir'),
        type='character',
        default=NULL,
        help='explicit out dir for imputation'
    ),
    make_option(
        c('--integrate'),
        type='logical',
        action='store_true',
        default=FALSE,
        help='whether or not to integrate'
    ),
    make_option(
        c('--plot_type'),
        type='character',
        default=NULL,
        help='what kind of plots to plot'
    ),
    make_option(
        c('--by_condition'),
        type='logical',
        action='store_true',
        default=FALSE,
        help='by condition or not'
    ),
    make_option(
        c('--species'),
        type='character',
        default=NULL,
        help='species'
    ),
    make_option(
        c('--mito_cutoff'),
        type='integer',
        default=NULL,
        help='upper umi cutoff'
    ),
    make_option(
        c('--upper_umi_cutoff'),
        type='integer',
        default=25000,
        help='upper umi cutoff'
    ),
    make_option(
        c('--rare_gene_cutoff'),
        type='integer',
        default=10,
        help='required number of cells expressing a gene to keep that gene'
    ),
    make_option(
        c('--min_uniq_gene_cutoff'),
        type='integer',
        default=350,
        help='required number of uniq genes for a cell'
    ),
    make_option(
        c('--width'),
        type='integer',
        default=10,
        help='plot width'
    ),
    make_option(
        c('--height'),
        type='integer',
        default=5,
        help='plot height'
    ),
#    make_option(
#        c('--xvar'),
#        type='character',
#        default=NULL,
#        help='xvar'
#    ),
#    make_option(
#        c('--yvar'),
#        type='character',
#        default=NULL,
#        help='yvar'
#    ),
    make_option(
        c('--xlab'),
        type='character',
        default=NULL,
        help='x axis label'
    ),
    make_option(
        c('--ylab'),
        type='character',
        default=NULL,
        help='y axis label'
    ),
    make_option(
        c('--facet_var'),
        type='character',
        default=NULL,
        help='subplots by factor'
    ),
    make_option(
        c('--resolution'),
        type='numeric',
        default=0.1,
        help='resolution ranges between 0 and 1. Higher values increase the number of clusters'
    ),
    make_option(
        c('--pcs'),
        type='numeric',
        default=50,
        help='number of principal components'
    ),
    make_option(
        c('--min_dist'),
        type='numeric',
        default=0.3,
        help='minimum distance for umap projections'
    ),
    make_option(
        c('--nn'),
        type='numeric',
        default=30,
        help='number of nearest neighbors for umap'
    ),
    make_option(
        c('--use_integrated'),
        type='logical',
        action='store_true',
        default=FALSE,
        help='wheter to use the integrated seurat obj instead of the imputed one'
    )
)


opt <- parse_args(OptionParser(option_list=option_list))
run_r <- opt$run_r
workflow <- opt$workflow
out_impute_dir <- opt$out_impute_dir
obj_path <- opt$obj_path
species <- opt$species
mito_cutoff <- opt$mito_cutoff
upper_umi_cutoff <- opt$upper_umi_cutoff
rare_gene_cutoff <- opt$rare_gene_cutoff
min_uniq_gene_cutoff <- opt$min_uniq_gene_cutoff
integrate <- opt$integrate
plot_type <- opt$plot_type
condition <- opt$by_condition
resolution <- opt$resolution
pcs <- opt$pcs
use_integrated <- opt$use_integrated
#xvar <- opt$xvar
#yvar <- opt$yvar
xlab <- opt$xlab
ylab <- opt$ylab
width <- opt$width
height <- opt$height
min_dist <- opt$min_dist
nn <- opt$nn



# error functions
check_species <- function (raw_args, mito_cutoff) {
    if (is.null(species)) {
        if (is.null(mito_cutoff)) {
            print('ERROR: no species provided (--species) AND no mito cutoff provided (--mito_cutoff)')
            quit(save='no')
        }
    } else {
        if (!(species %in% c('human', 'mouse')) & !(is.null(mito_cutoff))) {
            print('ERROR: no default mito cutoff for provided species (--mito_cutoff) and no mito cutoff provided')
            quit(save='no')
        }
        if (is.null(mito_cutoff)) {
            if (species == 'mouse') {
                raw_args <- paste0(raw_args, ' --mito_cutoff 10')
            } else {
                raw_args <- paste0(raw_args, ' --mito_cutoff 5')
            }
        }
    }
}


# workflows
run_preprocess <- function(mtx_dir, rare_gene_cutoff, mito_cutoff, upper_umi_cutoff) {
    sample_dirs <- fetch_mtx_dirs(mtx_dir)
    merged_data <- merge_all(sample_dirs)
    preprocessed_obj <- apply_rare_gene_filter(merged_data$obj, merged_data$sample_ids, rare_gene_cutoff)
    preprocessed_obj <- apply_min_uniq_gene_filter(preprocessed_obj, min_uniq_gene_cutoff)
    preprocessed_obj <- apply_upper_umi_cutoff(preprocessed_obj, upper_umi_cutoff)
    preprocessed_obj <- apply_mito_filter(preprocessed_obj, mito_cutoff)
    saveRDS(preprocessed_obj, file="./data/seurat-objects/preprocessed.rds")
    apply_imputation(preprocessed_obj)
    if (integrate) {
        apply_integration(preprocessed_obj)
    }
}


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


cluster <- function(out_dir, pcs, resolution, use_integrated) {
    if (use_integrated) {
        obj_path <- './data/seurat-objects/integrated.rds'
    } else {
        obj_path <- './data/seurat-objects/preprocessed.rds'
    }
    obj <- load_seurat_obj(obj_path)
    obj <- apply_dim_reduce(obj, pcs)
    apply_clustering(obj, out_dir, pcs, resolution, use_integrated)
}


### jobscripts ###
write_align_jobs <- function(raw_args) {
script <- c(
    "#!/bin/bash",
    "",
    "#SBATCH -A syb111",
    "#SBATCH -N $num_paired_files",
    "#SBATCH -t 6:00:00",
    "#SBATCH -J fastqc",
    "#SBATCH -o ./scripts/alignment/logs/fastqc.%J.out",
    "#SBATCH -e ./scripts/alignment/logs/fastqc.%J.err",
    "",
    "module load python",
    "echo fastqc_bin: $fastqc_bin",
    "echo",
    "srun -n $num_files python ./scripts/alignment/mpi_fastqc.py $fastqc_bin $data_dir $out_fastqc_dir",
    "",
    "source activate /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/python_andes",
    "srun -N 1 -n 1 multiqc $out_fastqc_dir -n fastqc_report.html -o $out_qc_dir --no-data-dir"
    )
    out_paths <- c("./scripts/jobs/fastqc.sbatch")
    write.table(script, file=out_path, sep="\n", row.names=FALSE, col.names=FALSE, quote=FALSE)

script <- c(
    "#!/bin/bash",
    "",
    "#SBATCH -A SYB111",
    paste0("#SBATCH -N ", num_paired_files),
    "#SBATCH -t 6:00:00",
    "#SBATCH -J STAR_align",
    "#SBATCH -o ./scripts/alignment/logs/align.%J.out",
    "#SBATCH -e ./scripts/alignment/logs/align.%J.err",
    "",
    "module load python",
    paste0("echo star_bin: ", star_bin),
    "echo",
    paste0("srun -n ", num_paired_files, "python ./scripts/alignment/mpi_align.py ", raw_args),
    "",
    "source activate /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/python_andes",
    paste0("mv ", out_align_dir, "/*Log*.out", out_align_dir, "/logs"),
    paste0("srun -N 1 -n 1 multiqc ", out_align_dir, "/logs -n align_report.html -o ", out_qc_dir, " --no-data-dir"),
    paste0("mv ", out_align_dir, "/*Solo.out ", out_mtx_dir)
    )
    out_paths <- c(out_paths, "./scripts/jobs/align.sbatch")
    write.table(script, file=out_path, sep="\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    return(out_paths)
}


#TODO: calculate time needed
write_preprocess_job <- function(raw_args) {
script <- c(
    "#!/bin/bash",
    "",
    "#SBATCH -A SYB111",
    "#SBATCH -N 1",
    "#SBATCH -t 2:00:00",
    "#SBATCH -J preprocess",
    "#SBATCH -o ./scripts/preprocess/logs/preprocess.%J.out",
    "#SBATCH -e ./scripts/preprocess/logs/preprocess.%J.err",
    "",
    #"source activate /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/deseq2_andes",
    "source activate /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_seurat",
    "",
    paste0("srun -n 1 Rscript ./sc-flow.R ", raw_args)
    )
    out_path <- "./scripts/jobs/preprocess.sbatch"
    write.table(script, file=out_path, sep="\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    return(c(out_path))
}

#
#TODO: calculate time needed
write_impute_job <- function(raw_args) {
script <- c(
    "#!/bin/bash",
    "",
    "#SBATCH -A SYB111",
    "#SBATCH -N 1",
    "#SBATCH -t 2:00:00",
    "#SBATCH -J impute",
    "#SBATCH -p gpu",
    "#SBATCH -o ./scripts/preprocess/logs/impute.%J.out",
    "#SBATCH -e ./scripts/preprocess/logs/impute.%J.err",
    "",
    "source activate /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/deseq2_andes",
    "",
    paste0("srun -n 1 Rscript ./sc-flow.R ", raw_args)
    )
    out_path <- "./scripts/jobs/impute.sbatch"
    write.table(script, file=out_path, sep="\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    return(c(out_path))
}


#TODO: calculate time needed
write_plot_job <- function(raw_args) {
script <- c(
    "#!/bin/bash",
    "",
    "#SBATCH -A SYB111",
    "#SBATCH -N 1",
    "#SBATCH -t 1:00:00",
    paste0("#SBATCH -J plot_", plot_type),
    paste0("#SBATCH -o ./scripts/plots/logs/plot_", plot_type, ".%J.out"),
    paste0("#SBATCH -e ./scripts/plots/logs/plot_", plot_type, ".%J.err"),
    "",
    "source activate /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/deseq2_andes",
    "",
    paste0("srun -n 1 Rscript ./sc-flow.R ", raw_args)
    )
    out_path <- paste0("./scripts/jobs/plot_", plot_type, ".sbatch")
    write.table(script, file=out_path, sep="\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    return(c(out_path))
}

submit_job <- function(path) {
    command <- paste0('sbatch ', path)
    system(command)
}


### WORKFLOW ###
valid_workflow_list <- c('align', 'preprocess', 'plot', 'impute', 'cluster')
valid_plot_type_list <- c('qc', 'umap')


if (is.null(workflow) | !(workflow %in% valid_workflow_list)) {
    print(paste0('ERROR: missing or invalid workflow specified (', workflow, ').'))
    quit(save='no')
}


if (!run_r) {
    raw_args <- paste0(raw_args, ' --run_r')
    if (workflow == 'align') {
        raw_args <- check_species(raw_args, mito_cutoff)
        job_paths <- write_align_jobs(raw_args)
    } else if (workflow == 'preprocess') {
        raw_args <- check_species(raw_args, mito_cutoff)
        job_paths <- write_preprocess_job(raw_args)
    } else if (workflow == 'impute') {
        job_paths <- write_impute_job(raw_args)
    } else if (workflow == 'plot') {
        if (is.null(plot_type) | !(plot_type %in% valid_plot_type_list)) {
            print(paste0('ERROR: missing or invalid plot type specified (', plot_type, ').'))
            quit(save='no')
        }
        if (plot_type == 'qc') {
            raw_args <- check_species(raw_args, mito_cutoff)
        }
        job_paths <- write_plot_job(raw_args)
    }
    # submit the job from the command line
    for (path in job_paths) {
        submit_job(path)
    }
### when the job script runs sc-flow.R ###
} else {
    source('./scripts/seurat/utils.R')
    if (workflow == 'preprocess') {
        source('./scripts/preprocess/preprocess.R')
        source('./scripts/preprocess/alra.R')
        run_preprocess(mtx_dir, rare_gene_cutoff, mito_cutoff, upper_umi_cutoff)
    } else if (workflow == 'impute') {
        source('./scripts/preprocess/preprocess.R')
        source('./scripts/preprocess/alra.R')
        obj <- load_seurat_obj(obj_path)
        apply_imputation(obj, out_impute_dir)
    } else if (workflow == 'cluster') {
        source('./scripts/seurat/cluster.R')
        cluster(meta_dir, pcs, resolution, use_integrated)
    } else if (workflow == 'plot') {
        source('./scripts/plots/plots.R')
        if (plot_type == 'qc') {
            plot_qc(metadata, xvar, yvar, condition, mito_cutoff,
                    plot_dir, plot_name, width, height)
        } else if (plot_type == 'umap') {
            source('./scripts/seurat/umap.R')
            gen_umap(meta_dir, pcs, min_dist, nn, use_integrated)
            plot_umap(meta_dir, plot_dir, width, height, pcs,
                      resolution, use_integrated, min_dist, nn)
        }
    }
}
