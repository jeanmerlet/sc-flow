# Libraries
suppressPackageStartupMessages({
    library(optparse)
})

raw_args <- paste0(commandArgs(trailingOnly = TRUE), collapse=' ')

mtx_dir <- './data/count-matrices/'
plot_dir <- './plots/'
meta_dir <- './data/metadata/'


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
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

run_r <- opt$run_r
workflow <- opt$workflow
species <- opt$species
mito_cutoff <- opt$mito_cutoff
upper_umi_cutoff <- opt$upper_umi_cutoff
rare_gene_cutoff <- opt$rare_gene_cutoff
integrate <- opt$integrate
plot_type <- opt$plot_type
condition <- opt$by_condition
#xvar <- opt$xvar
#yvar <- opt$yvar
xlab <- opt$xlab
ylab <- opt$ylab
width <- opt$width
height <- opt$height


# workflows
run_preprocess <- function(mtx_dir, rare_gene_cutoff, mito_cutoff, upper_umi_cutoff) {
    sample_dirs <- fetch_mtx_dirs(mtx_dir)
    merged_data <- merge_all(sample_dirs)
    preprocessed_obj <- apply_rare_gene_filter(merged_data$obj, merged_data$sample_ids, rare_gene_cutoff)
    preprocessed_obj <- apply_mito_filter(preprocessed_obj, species, mito_cutoff)
    preprocessed_obj <- apply_upper_umi_cutoff(preprocessed_obj, upper_umi_cutoff)
    apply_imputation(preprocessed_obj)
    if (integrate) {
        apply_integration(preprocessed_obj)
    }
}
#source('./scripts/preprocess/preprocess.R')
#source('./scripts/preprocess/alra.R')
#run_preprocess(mtx_dir, 10, 0.10, 25000)
#quit(save='no')


plot <- function(plot_type, meta_path, xvar, yvar, xlab, ylab, facet_var, 
                 condition, mito_cutoff, plot_dir, plot_name, width, height) {
    
}



# jobscripts
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
    "#SBATCH -t 6:00:00",
    "#SBATCH -J preprocess",
    "#SBATCH -o ./scripts/preprocess/logs/preprocess.%J.out",
    "#SBATCH -e ./scripts/preprocess/logs/preprocess.%J.err",
    "",
    "source activate /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/deseq2_andes",
    "",
    paste0("srun -n 1 Rscript ./sc-flow.R ", raw_args)
    )
    out_path <- "./scripts/jobs/preprocess.sbatch"
    write.table(script, file=out_path, sep="\n", row.names=FALSE, col.names=FALSE, quote=FALSE)
    return(c(out_path))
}


submit_job <- function(path) {
    command <- paste0('sbatch ', path)
    system(command)
}


valid_workflow_list <- c('align', 'preprocess', 'plot')
valid_plot_type_list <- c('qc')


if (is.null(workflow) | !(workflow %in% valid_workflow_list)) {
    print(paste0('ERROR: missing or invalid workflow specified (', workflow, ').'))
    quit(save='no')
}


if (!run_r) {
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
    }
    raw_args <- paste0(raw_args, ' --run_r')
    if (workflow == 'align') {
        job_paths <- write_align_jobs(raw_args)
    } else if (workflow == 'preprocess') {
        job_paths <- write_preprocess_job(raw_args)
    } else if (workflow == 'plot') {
        # make plotting job script
    }
    # submit the job from the command line
    for (path in job_paths) {
        submit_job(path)
    }
} else {
    if (workflow == 'preprocess') {
        source('./scripts/preprocess/preprocess.R')
        source('./scripts/preprocess/alra.R')
        run_preprocess(mtx_dir, rare_gene_cutoff, mito_cutoff, upper_umi_cutoff)
    } else if (workflow == 'cluster') {
        print('clustering')
    } else if (workflow == 'plot') {
        if (is.null(plot_type) | !(plot_type %in% valid_plot_type_list)) {
            print(paste0('ERROR: missing or invalid plot type specified (', plot_type, ').'))
            quit(save='no')
        }
        source('./scripts/plots/plots.R')
        if (plot_type == 'qc') {
            meta_path <- paste0(meta_dir, 'cell_meta_imputed.tsv')
            metadata <- read.table(meta_path, sep='\t')
            xvar <- 'sample_ids'
            yvar <- 'Mito'
            violin_plot(metadata, xvar, yvar, xlab, ylab, condition, mito_cutoff,
                        plot_dir, paste0(plot_type, '.png'), width, height)
        }
    }
}
