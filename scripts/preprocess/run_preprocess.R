# Libraries
suppressPackageStartupMessages({
    library(optparse)
})


mtx_dir <- './data/count-matrices/'


source('./scripts/preprocess/preprocess.R')
source('./scripts/preprocess/alra.R')

# Argument parser
option_list <- list(
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

species <- opt$species
mito_cutoff <- opt$mito_cutoff
upper_umi_cutoff <- opt$upper_umi_cutoff
rare_gene_cutoff <- opt$rare_gene_cutoff
integrate <- opt$integrate

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


# run preprocessing workflow
sample_dirs <- fetch_mtx_dirs(mtx_dir)
merged_data <- merge_all(sample_dirs)
preprocessed_obj <- apply_rare_gene_filter(merged_data$obj, merged_data$sample_ids, rare_gene_cutoff)
preprocessed_obj <- apply_mito_filter(preprocessed_obj, species, mito_cutoff)
preprocessed_obj <- apply_upper_umi_cutoff(preprocessed_obj, upper_umi_cutoff)
apply_imputation(preprocessed_obj)
if (integrate) {
    apply_integration(preprocessed_obj)
}
