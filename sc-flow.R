# Libraries
suppressPackageStartupMessages({
    library(optparse)
})

raw_args <- paste0(commandArgs(trailingOnly = TRUE), collapse=' ')

mtx_dir <- './data/count-matrices/'


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
        c('--species'),
        type='character',
        default=NULL,
        help='species'
    ),
    make_option(
        c('--integrate'),
        type='logical',
        action='store_true',
        default=FALSE,
        help='whether or not to integrate'
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


# jobscripts
write_alignment_jobs <- function() {

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
    write.table(script, file=out_path, sep = "\n",row.names = FALSE,col.names = FALSE, quote = FALSE)
    return(out_path)
}


submit_job <- function(path) {
    command <- paste0('sbatch ', path)
    system(command)
}


valid_workflow_list <- c('align', 'preprocess')

if (is.null(workflow) | !(workflow %in% valid_workflow_list)) {
    print('ERROR: no species provided (--species) AND no mito cutoff provided (--mito_cutoff)')
    print(paste0('ERROR: no or invalid job (', workflow, ') specified.'))
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
        job_path <- write_align_jobs()
    } else if (workflow == 'preprocess') {
        job_path <- write_preprocess_job(raw_args)
    }
    # run the job from the command line
    submit_job(job_path)
} else {
    if (workflow == 'preprocess') {
        source('./scripts/preprocess/preprocess.R')
        source('./scripts/preprocess/alra.R')
        run_preprocess(mtx_dir, rare_gene_cutoff, mito_cutoff, upper_umi_cutoff)
    } else if (workflow == 'cluster') {
    }
}
