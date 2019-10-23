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

# input feature gene symbol to ensembl mapping ----------------------------
feature.names = read.delim(file = "/diskmnt/Projects/ccRCC_scratch/snRNA_Processed_Data/inferCNV/inputs/features.tsv.gz", header = FALSE,stringsAsFactors = FALSE)

# input gene order file ----------------------------------------------------
tab_gene_order_ensembl = read.delim(file = path_gene_order_file, header = FALSE,stringsAsFactors = FALSE)

# input seurat object
seurat_object <- readRDS(file = "/diskmnt/Projects/ccRCC_scratch/snRNA_Processed_Data/analysis_results/integration/integrate_seurat_objects/20190927.v1/renal_integrated.20190927.v1.RDS")

# get cells 2 process -----------------------------------------------------
sample_bc <- rownames(seurat_object@meta.data)[seurat_object@meta.data$orig.ident == sample_id]
sample_bc

# get barcode 2 cluster table ---------------------------------------------
anno_tab <- seurat_object@meta.data[seurat_object@meta.data$orig.ident == sample_id,]
anno_tab$barcode <- rownames(anno_tab)
anno_tab <- anno_tab[, c("barcode", "seurat_clusters")]
path_anno_file <- paste0(dir_out, "annotation_file.", run_id, ".txt")
write.table(x = anno_tab, file = path_anno_file, quote = F, sep = "\t", row.names = F, col.names = F)

# get ref group names -----------------------------------------------------
all_group_names <- unique(anno_tab$seurat_clusters)
all_group_names
tumor_group_names <- c("1", "3", "4", "6", "7", "8")
ref_group_names <- all_group_names[!(all_group_names %in% tumor_group_names)]

if (!file.exists(paste0(dir_out, "infercnv_object.", run_id, ".RDS")) {
# get raw counts matrix
raw_exp_mat <- seurat_object@assays$RNA@counts[, sample_bc]
rownames(raw_exp_mat) <- plyr::mapvalues(rownames(raw_exp_mat), from = feature.names$V2,to = feature.names$V1)

# expression data - remove gene not mapped in the gene position file
missing_sym  <- rownames(raw_exp_mat)[!(rownames(raw_exp_mat) %in% tab_gene_order_ensembl$V1)]
## only missing ~750 genes
missing_row  <- match(missing_sym, rownames(raw_exp_mat))
raw_exp_mat_clean <- raw_exp_mat[-missing_row,]
raw_exp_mat_clean <- as.matrix(raw_exp_mat_clean)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_exp_mat_clean,
                                    annotations_file=path_anno_file,
                                    delim="\t",
                                    gene_order_file=path_gene_order_file,
                                    ref_group_names= ref_group_names)
saveRDS(object = infercnv_obj, file = paste0(dir_out, "infercnv_object.", run_id, ".RDS"))
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
