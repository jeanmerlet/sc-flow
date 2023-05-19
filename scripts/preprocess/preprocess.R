# Libraries
suppressPackageStartupMessages({
    library(Seurat)
    library(readr)
    library(stringr)
    library(DropletUtils)
    library(Matrix)
    library(R.utils)
    library(stringr)
})

code_dir <- './scripts/preprocess'
source(paste0(code_dir, '/alra.R'))

root_dir <- '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/sc-flow/data/count-matrices/'

#TODO: conditional Gene / GeneFull for introns
#TODO: conditional filtered / raw
# get a list of mtx dirs
fetch_mtx_dirs <- function(root_dir) {
    sample_dirs <- list.files(root_dir)
    sample_dirs <- do.call("c",lapply(sample_dirs, function(sample_dir){
        current_dir <- paste0(root_dir, sample_dir, "/Gene/", "filtered")
        f <- list.files(current_dir)
        lapply(c("barcodes.tsv","features.tsv","matrix.mtx"),function(file){
            if(!(paste0(file, ".gz") %in% f)) {
                gzip(paste0(current_dir, "/", file),paste0(current_dir, "/", file,".gz"),
                overwrite = TRUE,remove = FALSE)
            }})
        return(current_dir)
    }))
    return(sample_dirs)
}


#TODO: minutes vs seconds
# read all samples and combine into one
merge_all <- function(sample_dirs) {
    start_time <- Sys.time()
    first <- TRUE
    sample_ids <- c()
    for (sample_dir in sample_dirs) {
        mtx <- Read10X(sample_dir)
        sample_id <- str_match(sample_dir, "count-matrices/(.*)_Solo.out")[1,2]
        colnames(mtx) <- paste0(colnames(mtx), '_', sample_id)
        rownames(mtx) <- toupper(rownames(mtx))
        seurat_obj <- CreateSeuratObject(counts=mtx)
        sample_ids <- c(sample_ids, rep(sample_id, length(colnames(mtx))))
        if (first == TRUE) {
            combined_seurat_obj <- seurat_obj
            first <- FALSE
        } else {
            combined_seurat_obj <- merge(combined_seurat_obj, seurat_obj)
        }
    }
	print(paste0("Sample combining runtime: ",round(end_time - start_time,2)," minutes"))
    saveRDS(combined_seurat_obj, file="./data/seurat-objects/filtered.rds")
    return(list(combined_seurat_obj, sample_ids))
}


#TODO: deal with intronic
#TODO: implement manual mito cutoff setting
# filter cells by mitochondrial gene expression
apply_mito_filter <- function(obj, species, cutoff=NULL) {
    if (is.null(cutoff)) {
        if (species == 'mouse') {
            cutoff <- 10
        } else if (species == 'human') {
            cutoff <- 5
        } else {
            cutoff <- 10
            print(paste0('WARNING: mitochondrial cutoff for ', species,
                        ' unknown. Set to 10% by default. Set a cutoff with -insert argument-'))
        }
    }
    obj[["Mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj <- subset(obj, subset=(Mito < cutoff))
    return(obj)
}


#TODO: manual cutoff
apply_rare_gene_filter <- function(obj, sample_ids, cutoff=10) {
    counts <- GetAssayData(object=obj, slot="counts")
    keep_genes <- Matrix::rowSums(counts) > cutoff
    filtered_counts <- counts[keep_genes, ]
    obj <- CreateSeuratObject(filtered_counts)
    #obj <- AddMetaData(object=obj, metadata=sample_ids, col.name='sample_ids')
    obj@meta.data$sample_ids <- sample_ids
    test_list <- SplitObject(obj, split.by='sample_ids')
    return(obj)
}


#TODO: manual cutoff
apply_upper_umi_cutoff <- function(obj, cutoff=25000) {
    obj <- subset(obj, subset=(nCount_RNA < cutoff))
    return(obj)
}


#TODO: minutes vs seconds
apply_integration <- function(obj) {
	start_time <- Sys.time()
    obj_list <- SplitObject(obj, split.by='sample_ids')
    obj_list <- lapply(X=obj_list, FUN=function(x) {
        x <- NormalizeData(x, normalization.method='LogNormalize', scale.factor=10000, verbose=FALSE)
        x <- FindVariableFeatures(x)
    })
    features <- SelectIntegrationFeatures(object.list=obj_list)
    anchors <- FindIntegrationAnchors(object.list=obj_list, anchor.features=features)
    obj <- IntegrateData(anchorset=anchors)
    DefaultAssay(obj) <- 'integrated'
    saveRDS(obj, file="./data/seurat-objects/integrated.rds")
	end_time <- Sys.time()
	print(paste0("Integration runtime: ",round(end_time - start_time,2)," minutes"))
}


#TODO: minutes vs seconds
apply_imputation <- function(obj) {
	start_time <- Sys.time()
    obj <- NormalizeData(obj, normalization.method='LogNormalize', scale.factor=10000, verbose=FALSE)
	alra_result <- alra(t(as.matrix(GetAssayData(object=obj, slot = "data"))))
	obj_alra <- t(alra_result[[3]])
	colnames(obj_alra) <- rownames(obj@meta.data)
	obj_alra <- Matrix(obj_alra,sparse = T)
	obj <- SetAssayData(object = obj,slot = "data",new.data = obj_alra)
	saveRDS(obj, file = "./data/seurat-objects/imputed.rds")	
	end_time <- Sys.time()
	print(paste0("Imputation runtime: ",round(end_time - start_time,2)," minutes"))
}


sample_dirs <- fetch_mtx_dirs(root_dir)
obj_data <- merge_all(sample_dirs)
obj <- obj_data[[1]]
sample_ids <- obj_data[[2]]
pp_obj <- apply_rare_gene_filter(obj, sample_ids)
pp_obj <- apply_mito_filter(pp_obj, 'mouse')
pp_obj <- apply_upper_umi_cutoff(pp_obj, 1000)

imputation_test <- apply_imputation(pp_obj)

#integration_test <- apply_integration(pp_obj)