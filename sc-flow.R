# Libraries
suppressWarnings(suppressPackageStartupMessages({
    library(optparse)
    library(future)
}))
# /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/deseq2_andes
# /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_seurat

raw_args <- paste0(commandArgs(trailingOnly = TRUE),collapse=' ')

plan('multicore',workers = 8)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = 100 * 1024^3)

raw_dir <- './data/raw/'
fastqc <- '/lustre/orion/syb111/proj-shared/Tools/frontier/FastQC/fastqc'
qc_dir <- './qc/'
fastqc_out_dir <- paste0(qc_dir,'/reads')
star <- '/lustre/orion/syb111/proj-shared/Tools/frontier/STAR-2.7.9a/bin/Linux_x86_64/STAR'
bam_dir <- './data/bam/'
mtx_dir <- './data/count-matrices/'
obj_dir <- './data/seurat-objects/'
plot_dir <- './plots/'
meta_dir <- './data/metadata/'
de_dir <- './data/de/'


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
        c('--color_by'),
        type='character',
        default='none',
        help='what to color plots by. valid: sample_ids, condition, clusters'
    ),
    make_option(
        c('--split_by'),
        type='character',
        default='none',
        help='what to split plots by. valid: sample_ids, condition, clusters'
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
        c('--resolution'),
        type='numeric',
        default=0.3,
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
        c('--p_value'),
        type='numeric',
        default=0.1,
        help='adjusted p-value threshold for differential expression'
    ),
    make_option(
        c('--top_n_genes'),
        type='numeric',
        default=5,
        help='number of genes per diff_type for volcano plot'
    ),
    make_option(
        c('--p_cutoff'),
        type='numeric',
        default=0.05,
        help='adjusted p-value cutoff for volcano plot'
    ),
    make_option(
        c('--log2fc'),
        type='numeric',
        default=0.5,
        help='minimum absolute log2 fold change for volcano plot'
    ),
    make_option(
        c('--use_integrated'),
        type='logical',
        action='store_true',
        default=FALSE,
        help='whether to use the integrated seurat obj instead of the imputed one'
    ),
    make_option(
        c('--load_meta'),
        type='logical',
        action='store_true',
        default=FALSE,
        help='whether to load sample metadata from the user meta file (sample_meta.tsv)'
    ),
    make_option(
        c('--diff_type'),
        type='character',
        default='cluster',
        help='type of differential expression comparison to run'
    ),
    make_option(
        c('--condition_group'),
        type='character',
        default=NULL,
        help='group level for the condition variable'
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
color_by <- opt$color_by
split_by <- opt$split_by
resolution <- opt$resolution
pcs <- opt$pcs
use_integrated <- opt$use_integrated
xlab <- opt$xlab
ylab <- opt$ylab
width <- opt$width
height <- opt$height
min_dist <- opt$min_dist
nn <- opt$nn
p_value <- opt$p_value
diff_type <- opt$diff_type
log2fc <- opt$log2fc
top_n_genes <- opt$top_n_genes
p_cutoff <- opt$p_cutoff
load_meta <- opt$load_meta


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


### jobscripts ###
write_align_jobs <- function(raw_args) {
num_files <- length(list.files(raw_dir))
num_paired_files <- num_files/2
script <- c(
    "#!/bin/bash",
    "",
    "#SBATCH -A syb111",
    paste0("#SBATCH -N ",num_paired_files), # take raw_dir and calculate
    "#SBATCH -t 6:00:00",
    "#SBATCH -J fastqc",
    "#SBATCH -o ./scripts/alignment/logs/fastqc.%J.out",
    "#SBATCH -e ./scripts/alignment/logs/fastqc.%J.err",
    "",
    "module load python",
    paste0("echo fastqc_bin: ",fastqc), # placeholder
    "echo",
    paste0("srun -n ",num_files," python ./scripts/alignment/mpi_fastqc.py ",fastqc," ",fastqc_out_dir," ",fastqc_out_dir), # placeholders
    "",
    "source /lustre/orion/syb111/proj-shared/Tools/frontier/load_anaconda.sh",
    "conda activate sc-flow",
    "",
    paste0("srun -N 1 -n 1 multiqc ",fastqc_out_dir," -n fastqc_report.html -o ",fastqc_out_dir," --no-data-dir") # add path to raw_dir, qc_dir, and raw_args
    )
    out_path <- "./scripts/jobs/fastqc.sbatch"
    out_paths <- c(out_path)
    write.table(script,file = out_path,sep = "\n",row.names = FALSE,col.names = FALSE,quote = FALSE)

script <- c(
    "#!/bin/bash",
    "",
    "#SBATCH -A SYB111",
    paste0("#SBATCH -N ",num_paired_files), # take raw_dir and calculate
    "#SBATCH -t 6:00:00",
    "#SBATCH -J STAR_align",
    "#SBATCH -o ./scripts/alignment/logs/align.%J.out",
    "#SBATCH -e ./scripts/alignment/logs/align.%J.err",
    "",
    "module load python",
    paste0("echo star_bin: ",star), # placeholder 
    "echo",
    paste0("srun -n ",num_paired_files,"python ./scripts/alignment/mpi_align.py ",raw_args), # placeholders
    "",
    "source /lustre/orion/syb111/proj-shared/Tools/frontier/load_anaconda.sh",
    "conda activate sc-flow",
    "",
    paste0("mv ",bam_dir,"/*Log*.out",bam_dir,"/logs"), # placeholders
    paste0("srun -N 1 -n 1 multiqc ",bam_dir,"/logs -n align_report.html -o ",qc_dir," --no-data-dir"), # placeholders
    paste0("mv ",bam_dir,"/*Solo.out ",mtx_dir) # placeholders
    )
    out_path <- "./scripts/jobs/align.sbatch"
    out_paths <- c(out_paths,out_path)
    write.table(script,file = out_path,sep = "\n",row.names = FALSE,col.names = FALSE,quote = FALSE)
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
    "source /lustre/orion/syb111/proj-shared/Tools/frontier/load_anaconda.sh",
    "conda activate sc-flow",
    "",
    paste0("srun -n 1 Rscript ./sc-flow.R ",raw_args)
    )
    out_path <- "./scripts/jobs/preprocess.sbatch"
    write.table(script,file = out_path,sep = "\n",row.names = FALSE,col.names = FALSE,quote = FALSE)
    return(c(out_path))
}


#TODO: calculate time needed
write_impute_job <- function(raw_args) {
script <- c(
    "#!/bin/bash",
    "",
    "#SBATCH -A SYB111",
    "#SBATCH -N 1",
    "#SBATCH -t 2:00:00",
    "#SBATCH -J impute",
    "#SBATCH -o ./scripts/preprocess/logs/impute.%J.out",
    "#SBATCH -e ./scripts/preprocess/logs/impute.%J.err",
    "",
    "source /lustre/orion/syb111/proj-shared/Tools/frontier/load_anaconda.sh",
    "conda activate sc-flow",
    "",
    paste0("srun -n 1 Rscript ./sc-flow.R ",raw_args)
    )
    out_path <- "./scripts/jobs/impute.sbatch"
    write.table(script,file = out_path,sep = "\n",row.names = FALSE,col.names = FALSE,quote = FALSE)
    return(c(out_path))
}


#TODO: calculate time needed
write_cluster_job <- function(raw_args) {
script <- c(
    "#!/bin/bash",
    "",
    "#SBATCH -A SYB111",
    "#SBATCH -N 1",
    "#SBATCH -t 2:00:00",
    "#SBATCH -J cluster",
    "#SBATCH -o ./scripts/seurat/logs/cluster.%J.out",
    "#SBATCH -e ./scripts/seurat/logs/cluster.%J.err",
    "",
    "source /lustre/orion/syb111/proj-shared/Tools/frontier/load_anaconda.sh",
    "conda activate sc-flow",
    "",
    paste0("srun -n 1 Rscript ./sc-flow.R ",raw_args)
    )
    out_path <- "./scripts/jobs/cluster.sbatch"
    write.table(script,file = out_path,sep = "\n",row.names = FALSE,col.names = FALSE,quote = FALSE)
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
    paste0("#SBATCH -J plot_",plot_type),
    paste0("#SBATCH -o ./scripts/plots/logs/plot_",plot_type,".%J.out"),
    paste0("#SBATCH -e ./scripts/plots/logs/plot_",plot_type,".%J.err"),
    "",
    "source /lustre/orion/syb111/proj-shared/Tools/frontier/load_anaconda.sh",
    "conda activate sc-flow",
    "",
    paste0("srun -n 1 Rscript ./sc-flow.R ",raw_args)
    )
    out_path <- paste0("./scripts/jobs/plot_",plot_type,".sbatch")
    write.table(script,file = out_path,sep = "\n",row.names = FALSE,col.names = FALSE,quote = FALSE)
    return(c(out_path))
}

submit_job <- function(path) {
    command <- paste0('sbatch ',path)
    system(command)
}


#TODO: calculate time needed
write_diff_exp_job <- function(raw_args) {
script <- c(
    "#!/bin/bash",
    "",
    "#SBATCH -A SYB111",
    "#SBATCH -N 1",
    "#SBATCH -t 4:00:00",
    paste0("#SBATCH -J diff-exp_",diff_type),
    paste0("#SBATCH -o ./scripts/seurat/logs/diff-exp_",diff_type,".%J.out"),
    paste0("#SBATCH -e ./scripts/seurat/logs/diff-exp_",diff_type,".%J.err"),
    "",
    "source /lustre/orion/syb111/proj-shared/Tools/frontier/load_anaconda.sh",
    "conda activate sc-flow",
    "",
    paste0("srun -n 1 -c 8 Rscript ./sc-flow.R ",raw_args)
    )
    out_path <- paste0("./scripts/jobs/diff-exp_",diff_type,".sbatch")
    write.table(script,file = out_path,sep = "\n",row.names = FALSE,col.names = FALSE,quote = FALSE)
    return(c(out_path))
}

submit_job <- function(path) {
    command <- paste0('sbatch ',path)
    system(command)
}


### WORKFLOW ###
valid_workflow_list <- c('align','preprocess','plot','impute','cluster','diff_exp')
valid_plot_type_list <- c('qc','umap','volcano')
valid_diff_type_list <- c('cluster','condition')


if (is.null(workflow) | !(workflow %in% valid_workflow_list)) {
    print(paste0('ERROR: missing or invalid workflow specified (',workflow,').'))
    quit(save = 'no')
}


if (!run_r) {
    raw_args <- paste0(raw_args, ' --run_r')
    if (workflow == 'align') {
        raw_args <- check_species(raw_args,mito_cutoff)
        job_paths <- write_align_jobs(raw_args)
    } else if (workflow == 'preprocess') {
        raw_args <- check_species(raw_args,mito_cutoff)
        job_paths <- write_preprocess_job(raw_args)
    } else if (workflow == 'impute') {
        job_paths <- write_impute_job(raw_args)
    } else if (workflow == 'cluster') {
        job_paths <- write_cluster_job(raw_args)
    } else if (workflow == 'plot') {
#TODO: move to error function and have work when doing interactively with --run_r
        if (is.null(plot_type) | !(plot_type %in% valid_plot_type_list)) {
            print(paste0('ERROR: missing or invalid plot type specified (',plot_type,').'))
            quit(save = 'no')
        }
        if (plot_type == 'qc') {
            raw_args <- check_species(raw_args,mito_cutoff)
        }
        job_paths <- write_plot_job(raw_args)
    } else if (workflow == 'diff_exp') {
#TODO: do a check for diff_type and move that check to error function
	if(is.null(diff_type) | !(diff_type %in% valid_diff_type_list)) {
	    print(paste0('ERROR: missing or invalid diff type specified (',diff_type,').'))
            quit(save = 'no')
	}
        job_paths <- write_diff_exp_job(raw_args)
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
        run_preprocess(mtx_dir,meta_dir,rare_gene_cutoff,mito_cutoff,upper_umi_cutoff,load_meta)
    } else if (workflow == 'impute') {
        source('./scripts/preprocess/preprocess.R')
        source('./scripts/preprocess/alra.R')
        obj <- load_seurat_obj(obj_path)
        apply_imputation(obj, out_impute_dir)
    } else if (workflow == 'cluster') {
        source('./scripts/seurat/cluster.R')
        cluster(meta_dir,pcs,resolution,use_integrated)
    } else if (workflow == 'plot') {
        source('./scripts/plots/plots.R')
        if (plot_type == 'qc') {
            plot_qc(meta_dir,mito_cutoff,plot_dir,width,height)
        } else if (plot_type == 'umap') {
            source('./scripts/seurat/umap.R')
            gen_umap(meta_dir,pcs,min_dist,nn,use_integrated)
            plot_umap(meta_dir,plot_dir,width,height,pcs,resolution,use_integrated,min_dist,nn,color_by,split_by)
        } else if (plot_type == 'volcano') {
            plot_volcano(de_dir,plot_dir,diff_type,p_value,p_cutoff,log2fc,top_n_genes,width,height)
        }
    } else if (workflow == 'diff_exp') {
        source('./scripts/seurat/diff_exp.R')
        obj_path <- './data/seurat-objects/imputed.rds'
        obj <- load_seurat_obj(obj_path)
        #TODO: move integrated / preprocessed logic to utils
        if (use_integrated) {
            obj_type <- 'integrated'
        } else {
            obj_type <- 'preprocessed'
        }
        meta_path <- paste0(meta_dir,'clusters_obj-type-',obj_type,'_resolution-',resolution,'_num-pcs-',pcs,'.tsv')
        obj <- add_metadata(obj,meta_path)
        run_diff_exp(obj,de_dir,diff_type,condition_group,p_value)
    }
}
