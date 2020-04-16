#!/bin/sh

# Remove old logs
rm -r *html;
rm .nextflow.log*;

export REPO_MODULE="${HOME}/repo/sc_nextflow/pipelines/0025-qc_cluster"
export STUDY_DIR="${HOME}/repo/sc_nextflow/studies/gut-freeze002/ti-cd_healthy"

# Nextflow settings
export NXF_HOME=$(pwd)
export NXF_WORK="${NXF_HOME}/.nexflow_work"
export NXF_TEMP="${NXF_HOME}/.nexflow_temp"
export NXF_CONDA_CACHEDIR="${NXF_HOME}/.nexflow_conda"
export NXF_SINGULARITY_CACHEDIR="${NXF_HOME}/.nexflow_singularity"

# Farm specific settings
export QT_QPA_PLATFORM='offscreen'
export LSB_DEFAULT_JOBGROUP="/${USER}/nf"
export LSB_DEFAULTGROUP="team152"

nextflow run "${REPO_MODULE}/main.nf" \
         -profile "lsf" \
         --file_paths_10x "${STUDY_DIR}/file_paths_10x.tsv" \
         --file_metadata "${HOME}/repo/scrna_cellranger/sync_status/samples_metainfo.tsv" \
         --file_sample_qc "${STUDY_DIR}/params-sample_qc-no_filters.yml" \
         --genes_exclude_hvg "${STUDY_DIR}/../filters-variable_gene.tsv" \
         --genes_score "${REPO_MODULE}/example_runtime_setup/genes_score_v001.tsv" \
         --output_dir "nf-sample_qc=no_filters-parameter_sweep_v001" \
         -params-file "${STUDY_DIR}/params-analysis-parameter_sweep_v001.yml" \
         --run_multiplet \
         -resume
