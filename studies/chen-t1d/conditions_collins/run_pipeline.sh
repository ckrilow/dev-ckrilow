#!/bin/sh

# Remove old logs
# rm -r *html;
# rm .nextflow.log*;

export REPO_MODULE="/cluster/ifs/projects/collins/taylorhj/projects/sc_nextflow/pipelines/0025-qc_cluster"
export STUDY_DIR="/cluster/ifs/projects/collins/taylorhj/projects/sc_nextflow/studies/chen-t1d/conditions_collins"

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

# Shouldn't need -- from Anderson LSF configuration
#export LSB_DEFAULT_JOBGROUP="/${USER}/nf"
#export LSB_DEFAULTGROUP="team152"

# One may need to unset R_LIBS to get the conda install
# unset R_LIBS
# unset R_LIBS_USER

/home/taylorhj/miniconda3/envs/nextflow20/bin/nextflow run \
    "${REPO_MODULE}/main.nf" \
     -profile "sge" \
     --file_paths_10x "${STUDY_DIR}/file_paths_10x.tsv" \
     --file_metadata "${STUDY_DIR}/conditions_clinical_cr_metrics.tsv" \
     --file_sample_qc "${STUDY_DIR}/params-sample_qc-no_filters.yml" \
     --genes_exclude_hvg "${STUDY_DIR}/../data-variable_gene_filter_ribo_mt.tsv" \
     --genes_score "${STUDY_DIR}/../data-gene_scores.tsv" \
     --output_dir "${STUDY_DIR}/results/nf-sample_qc=no_filter_singlet-parameter_sweep_v001" \
     -params-file "${STUDY_DIR}/params-analysis-parameter_sweep_v001.yml" \
     --run_multiplet \
     -with-report \
     -with-trace \
     -with-timeline \
     -with-dag flowchart.png \
     -resume