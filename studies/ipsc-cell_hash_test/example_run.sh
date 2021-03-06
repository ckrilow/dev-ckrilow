#!/bin/sh

# Remove old logs
# rm -r *html;
# rm .nextflow.log*;

export REPO_MODULE="${HOME}/repo/sc_nextflow/pipelines/0025-qc_cluster"
export STUDY_DIR="${HOME}/repo/sc_nextflow/studies/ipsc-cell_hash_test"

# Nextflow settings
export NXF_OPTS="-Xms10G -Xmx10G"
# Uncomment this if get strange bus errors
# export NXF_OPTS="${NXF_OPTS} -Dleveldb.mmap=false" # No resume functionality
export NXF_HOME=$(pwd)
export NXF_WORK="${NXF_HOME}/.nextflow_work"
export NXF_TEMP="${NXF_HOME}/.nextflow_temp"
export NXF_CONDA_CACHEDIR="${NXF_HOME}/.nextflow_conda"
export NXF_SINGULARITY_CACHEDIR="${NXF_HOME}/.nextflow_singularity"

# Farm specific settings
export QT_QPA_PLATFORM='offscreen'
export LSB_DEFAULT_JOBGROUP="/${USER}/nf"
export LSB_DEFAULTGROUP="team152"

/software/hgi/installs/anaconda3/envs/nextflow20/bin/nextflow run \
    "${REPO_MODULE}/main.nf" \
     -profile "lsf" \
     --file_paths_10x "${STUDY_DIR}/file_paths_10x.tsv" \
     --file_metadata "${STUDY_DIR}/file_metadata.tsv" \
     --file_sample_qc "${STUDY_DIR}/params-sample_qc-no_filters.yml" \
     --genes_exclude_hvg "${STUDY_DIR}/filters-variable_gene.tsv" \
     --genes_score "${REPO_MODULE}/example_runtime_setup/genes_score_v001.tsv" \
     --output_dir "$(pwd)/nf-sample_qc=none-parameter_sweep_v001" \
     -params-file "${STUDY_DIR}/params-analysis-parameter_sweep_v001.yml" \
     --run_multiplet \
     -with-report \
     -with-trace \
     -with-timeline \
     -with-dag flowchart.png \
     -resume
