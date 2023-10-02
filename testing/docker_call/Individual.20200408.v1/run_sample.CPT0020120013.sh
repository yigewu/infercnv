#!/bin/bash

snRNA_aliquot_id=CPT0020120013
analysis_mode=subclusters
dir_run=/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/
docker_image_name=singlecellportal/infercnv
dir_infercnv=${dir_run}Resources/snRNA_Processed_Data/InferCNV/
dir_infercnv_outputs=${dir_infercnv}outputs/
dir_infercnv_inputs=${dir_infercnv}inputs/
run_id=20200408.v1
dir_infercnv_outputs_by_run=${dir_infercnv_outputs}Individual.${run_id}/
mkdir -p ${dir_infercnv_outputs_by_run}
dir_output=${dir_infercnv_outputs_by_run}${snRNA_aliquot_id}/
mkdir -p ${dir_output}
dir_raw_count_file=${dir_infercnv_inputs}raw_counts_matrix/
path_raw_count_file=${dir_raw_count_file}${snRNA_aliquot_id}.RNA_Count.tsv
dir_annotation_file=${dir_infercnv_inputs}annotations_file/Individual.${run_id}/
path_annotation_file=${dir_annotation_file}${snRNA_aliquot_id}.Barcode_Annotation.txt
path_log_file=${dir_output}${snRNA_aliquot_id}.$(date +%Y%m%d%H%M%S).log
path_gene_order_file=${dir_infercnv_inputs}gencode_v21_gen_pos.MatchCellRangerFeatures.NoDuplicates.20191005.v1.txt
cutoff=0.04
ref_group_names=Ref
num_ref_groups=20
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
        --num_ref_groups=${num_ref_groups} \
        --ref_group_names=${ref_group_names} &> ${path_log_file}&





