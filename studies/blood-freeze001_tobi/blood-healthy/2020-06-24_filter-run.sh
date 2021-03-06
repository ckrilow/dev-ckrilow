#!/bin/sh

# Remove old logs
# rm -r *html;
# rm .nextflow.log*;

# Custom to run
export REPO_MODULE="${PROJECT_HOME}/sc_nextflow/pipelines/0025-qc_cluster"
export STUDY_DIR="${PROJECT_HOME}/sc_nextflow/studies/blood-freeze001_tobi/blood-healthy"
export QC_FILENAME="params-sample_qc-mito0pt80_ngene150_singlet.yml"
export PARAMS_FILENAME="params-analysis-blood-custom.yml"

# Fixed
export OUTPUT_DIR="$(pwd)/nf-${PARAMS_FILENAME%.*}${QC_FILENAME%.*}"

# Nextflow settings
export NXF_OPTS="-Xms25G -Xmx25G"
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

# One may need to unset R_LIBS to get the conda install
# unset R_LIBS
# unset R_LIBS_USER

/software/hgi/installs/anaconda3/envs/nextflow20/bin/nextflow run \
    "${REPO_MODULE}/main.nf" \
     -profile "lsf" \
     --file_paths_10x "${STUDY_DIR}/file_paths_10x.tsv" \
     --file_metadata "${PROJECT_HOME}/scrna_cellranger/sync_status/samples_metainfo.tsv" \
     --file_sample_qc "${STUDY_DIR}/${QC_FILENAME}" \
     --genes_exclude_hvg "${STUDY_DIR}/../filters-variable_gene.tsv" \
     --genes_score "${REPO_MODULE}/example_runtime_setup/genes_score_v001.tsv" \
     --output_dir "${OUTPUT_DIR}" \
     -params-file "${STUDY_DIR}/${PARAMS_FILENAME}" \
     --run_multiplet \
     -with-report \
     -with-trace \
     -with-timeline \
     -with-dag flowchart.png \
     -resume
