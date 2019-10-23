#!/usr/bin/env Rscript

library("infercnv")

# set parameters ----------------------------------------------------------

version_tmp <- 1
sample_ids <- c("CPT0075130004_notFACS", "CPT0086820004_notFACS", "CPT0075140002", "CPT0001260013", "CPT0086350004")
sample_id <- "CPT0086350004"
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out_all <- "/diskmnt/Projects/ccRCC_scratch/snRNA_Processed_Data/inferCNV/outputs/"
dir_out <- paste0(dir_out_all, sample_id, "/")
dir.create(dir_out)
path_gene_order_file <- "/diskmnt/Projects/ccRCC_scratch/snRNA_Processed_Data/inferCNV/inputs/gencode_v21_gen_pos.ensembl_gene_id.txt" 

# input gene order file ----------------------------------------------------
tab_gene_order_ensembl = read.delim(file = path_gene_order_file, header = FALSE,stringsAsFactors = FALSE)

# path to annotation file
path_anno_file <- "/diskmnt/Projects/ccRCC_scratch/analysis_results/copy_number/run_infercnv/20191002.v1/annotation_file.20191002.v1.txt"

# get ref group names -----------------------------------------------------
ref_group_names <- c("C0", "C11", "C2", "C14", "C5", "C9", "C13", "C15", "C10", "C16", "C12")

if (!file.exists(paste0(dir_out, "infercnv_object.", run_id, ".RDS"))) {
# get matrix directory
dir_cellranger_output <- "/diskmnt/Projects/ccRCC_scratch/snRNA_Processed_Data/Cell_Ranger/outputs/"
dir_raw <- paste0(dir_cellranger_output, sample_id, "/outs/raw_feature_bc_matrix/")

# get direct paths to data
path_matrix <- paste0(dir_raw, "matrix.mtx.gz")
path_barcodes <- paste0(dir_raw, "barcodes.tsv.gz")
path_features <- paste0(dir_raw, "features.tsv.gz")

# read in matrix
raw_exp_mat <- readMM(file = path_matrix)
feature.names = read.delim(path_features, header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(path_barcodes, header = FALSE,stringsAsFactors = FALSE)

colnames(raw_exp_mat) = barcode.names$V1
rownames(raw_exp_mat) = feature.names$V1

# expression data - remove gene not mapped in the gene position file
missing_sym  <- rownames(raw_exp_mat)[!(rownames(raw_exp_mat) %in% tab_gene_order_ensembl$V1)]
missing_sym
## only missing ~750 genes, romove them
clean_exp_mat <- raw_exp_mat[!(rownames(raw_exp_mat) %in% missing_sym),]
rm(raw_exp_mat)
clean_exp_mat <- as.matrix(clean_exp_mat)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_exp_mat_clean,
                                    annotations_file=path_anno_file,
                                    delim="\t",
                                    gene_order_file=path_gene_order_file,
                                    ref_group_names= ref_group_names)
saveRDS(object = infercnv_obj, file = paste0(dir_out, "infercnv_object.", run_id, ".RDS"))
} else {
	infercnv_obj <- readRDS(paste0(dir_out, "infercnv_object.", run_id, ".RDS"))
}

for (i in c(0.1)) {	
	dir_out_i=paste0(dir_out, "raw_counts.cutoff.",i)
# perform infercnv operations to reveal cnv signal
	infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=i, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=dir_out_i, 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             mask_nonDE_genes = T,
                             include.spike=T  # used for final scaling to fit range (0,2) centered at 1.
                             )
}
