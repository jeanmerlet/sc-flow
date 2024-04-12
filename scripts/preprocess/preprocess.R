# Libraries
suppressPackageStartupMessages({
    library(Seurat)
    library(readr)
    library(stringr)
    library(DropletUtils)
    library(Matrix)
    library(R.utils)
})


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
merge_all <- function(sample_dirs, meta_dir, load_meta) {
    start_time <- Sys.time()
    if (load_meta) {
        user_meta <- read.table(paste0(meta_dir, 'sample_meta.tsv'), sep='\t', header=TRUE)
    }
    first <- TRUE
    sample_ids <- c()
    for (sample_dir in sample_dirs) {
        mtx <- Read10X(sample_dir)
        sample_id <- str_match(sample_dir, "count-matrices/(.*)_Solo.out")[1,2]
        if (load_meta) {
            condition <- user_meta[user_meta$sample_name==sample_id, ]$condition
        }
        colnames(mtx) <- paste0(colnames(mtx), '_', sample_id)
        rownames(mtx) <- toupper(rownames(mtx))
        seurat_obj <- CreateSeuratObject(counts=mtx)
        sample_ids <- c(sample_ids, rep(sample_id, length(colnames(mtx))))
        if (first == TRUE) {
            combined_seurat_obj <- seurat_obj
            if (load_meta) {
                row_names <- Cells(seurat_obj)
                obj_meta <- data.frame('condition'=rep(condition, length(row_names)))
                rownames(obj_meta) <- row_names
            }
            first <- FALSE
        } else {
            combined_seurat_obj <- merge(combined_seurat_obj, seurat_obj)
            if (load_meta) {
                row_names <- c(rownames(obj_meta), Cells(seurat_obj))
                obj_meta <- rbind(obj_meta, data.frame('condition'=rep(condition, length(Cells(seurat_obj)))))
                rownames(obj_meta) <- row_names
            }
        }
    }
	end_time <- Sys.time()
	print(paste0("Sample combining runtime: ",round(end_time - start_time, 2)," minutes"))
    saveRDS(combined_seurat_obj, file="./data/seurat-objects/filtered.rds")
    if (load_meta) {
        merged_data <- list(combined_seurat_obj, sample_ids, obj_meta)
        names(merged_data) <- c('obj', 'sample_ids', 'condition')
    } else {
        merged_data <- list(combined_seurat_obj, sample_ids)
        names(merged_data) <- c('obj', 'sample_ids')
    }
    return(merged_data)
}


#TODO: manual cutoff
apply_rare_gene_filter <- function(obj, sample_ids, cutoff, condition, load_meta) {
    counts <- GetAssayData(object=obj, slot="counts")
    keep_genes <- Matrix::rowSums(counts) > cutoff
    filtered_counts <- counts[keep_genes, ]
    obj <- CreateSeuratObject(filtered_counts)
    if (load_meta) {
        obj <- AddMetaData(obj, condition, col.name='condition')
    }
    obj <- AddMetaData(obj, sample_ids, col.name='sample_ids')
    #obj@meta.data$sample_ids <- sample_ids
    return(obj)
}


#TODO: deal with intronic
#TODO: implement manual mito cutoff setting
# filter cells by mitochondrial gene expression
apply_mito_filter <- function(obj, cutoff) {
    obj[["Mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj <- subset(obj, subset=(Mito < cutoff))
    write.table(obj@meta.data, file="./data/metadata/cell_meta_filtered.tsv", sep="\t")
    write.table(obj$RNA@meta.features, file="./data/metadata/gene_meta_filtered.tsv", sep="\t")
    return(obj)
}


apply_min_uniq_gene_filter <- function(obj, cutoff) {
    obj <- subset(obj, subset=(nFeature_RNA > cutoff))
    return(obj)
}


apply_upper_umi_cutoff <- function(obj, cutoff) {
    obj <- subset(obj, subset=(nCount_RNA < cutoff))
    return(obj)
}


#TODO: minutes vs seconds
apply_imputation <- function(obj, out_dir=NULL) {
	start_time <- Sys.time()
    obj <- NormalizeData(obj, normalization.method='LogNormalize', scale.factor=10000, verbose=FALSE)
	alra_result <- alra(t(as.matrix(GetAssayData(object=obj, slot = "data"))))
	obj_alra <- t(alra_result[[3]])
	colnames(obj_alra) <- rownames(obj@meta.data)
	obj_alra <- Matrix(obj_alra,sparse = T)
	obj <- SetAssayData(object = obj,slot = "data",new.data = obj_alra)
    if (is.null(out_dir)) {
        saveRDS(obj, file="./data/seurat-objects/imputed.rds")	
        write.table(obj@meta.data, file="./data/metadata/cell_meta_imputed.tsv", sep="\t")
        write.table(obj$RNA@meta.features, file="./data/metadata/gene_meta_imputed.tsv", sep="\t")
    } else {
        saveRDS(obj, file=paste0(out_dir, "imputed.rds"))
        write.table(obj@meta.data, file=paste0(out_dir, "cell_meta_imputed.tsv"), sep="\t")
        write.table(obj$RNA@meta.features, file=paste0(out_dir, "gene_meta_imputed.tsv"), sep="\t")
    }
	end_time <- Sys.time()
	print(paste0("Imputation runtime: ",round(end_time - start_time, 2)," minutes"))
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
    print(paste0("Integration runtime: ",round(end_time - start_time, 2)," minutes"))
}


# preprocess workflow
run_preprocess <- function(mtx_dir, meta_dir, rare_gene_cutoff, mito_cutoff, upper_umi_cutoff, load_meta) {
    sample_dirs <- fetch_mtx_dirs(mtx_dir)
    merged_data <- merge_all(sample_dirs, meta_dir, load_meta)
    obj <- merged_data$obj
    sample_ids <- merged_data$sample_ids
    if (load_meta) {
        condition <- merged_data$condition
    } else {
        condition <- NULL
    }
    preprocessed_obj <- apply_rare_gene_filter(obj, sample_ids, rare_gene_cutoff, condition, load_meta)
    preprocessed_obj <- apply_min_uniq_gene_filter(preprocessed_obj, min_uniq_gene_cutoff)
    preprocessed_obj <- apply_upper_umi_cutoff(preprocessed_obj, upper_umi_cutoff)
    preprocessed_obj <- apply_mito_filter(preprocessed_obj, mito_cutoff)
    saveRDS(preprocessed_obj, file="./data/seurat-objects/preprocessed.rds")
    apply_imputation(preprocessed_obj)
    if (integrate) {
        apply_integration(preprocessed_obj)
    }
}
