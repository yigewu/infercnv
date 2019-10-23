#!/bin/bash

analysis_mode=subclusters
dir_run=/diskmnt/Projects/ccRCC_scratch/
docker_image_name=singlecellportal/infercnv
integration_id=20191015.v1
dir_infercnv_input=/diskmnt/Projects/ccRCC_scratch/Resources/snRNA_Processed_Data/InferCNV/outputs/integration.${integration_id}/
mkdir -p ${dir_infercnv_input}

snRNA_aliquot_id=CPT0086350004
dir_output=${dir_infercnv_input}${snRNA_aliquot_id}_${analysis_mode}/
mkdir -p ${dir_output}

dir_raw_count_file=${dir_infercnv_input}
path_raw_count_file=${dir_raw_count_file}${snRNA_aliquot_id}/${snRNA_aliquot_id}.RNA_Count.20191016.v1.tsv
dir_annotation_file=${dir_raw_count_file}
path_annotation_file=${dir_annotation_file}${snRNA_aliquot_id}/${snRNA_aliquot_id}.Barcode_Annotation.20191016.v1.txt
path_log_file=${dir_output}${snRNA_aliquot_id}.$(date +%Y%m%d%H%M%S).log
path_gene_order_file=/diskmnt/Projects/ccRCC_scratch/Resources/snRNA_Processed_Data/InferCNV/inputs/gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.20191005.v1.txt
cutoff=0.04
ref_group_names="3,4,5,8,9,10,12,13,14,15,16"

docker run -v ${dir_run}:${dir_run} ${docker_image_name} inferCNV.R \
	--analysis_mode=${analysis_mode} \
	--raw_counts_matrix=${path_raw_count_file} \
	--annotations_file=${path_annotation_file} \
	--gene_order_file=${path_gene_order_file} \
	--cutoff=${cutoff} \
	--out_dir=${dir_output} \
	--cluster_by_groups \
	--denoise \
	--HMM \
	--num_threads=20 \
	--ref_group_names=${ref_group_names} &> ${path_log_file}&
