#!/bin/sh

source activate sc_preprocessing

# Remove old logs but not the most previous run
rm -r *html.*;
rm .nextflow.log.*;
rm flowchart.png.*;
rm trace.txt.*;

export REPO_MODULE="/cluster/ifs/projects/collins/taylorhj/projects/sc_nextflow/pipelines/0015-preprocessing"
export STUDY_DIR="/cluster/ifs/projects/collins/taylorhj/projects/sc_nextflow/studies/chen-t1d/conditions_collins"
#export OUTPUT_DIR="${STUDY_DIR}/results/min_parameter_sweep"
export OUTPUT_DIR="${STUDY_DIR}/results/nf-preprocessing=parker_singlet"

# Nextflow settings
export JAVA_HOME="/cluster/ifs/projects/collins/taylorhj/bin/jre1.8.0_251"
export JAVA_CMD="/cluster/ifs/projects/collins/taylorhj/bin/jre1.8.0_251/bin/java"
export NXF_OPTS="-Xms25G -Xmx25G"
# Uncomment this if get strange bus errors
# export NXF_OPTS="${NXF_OPTS} -Dleveldb.mmap=false" # No resume functionality
export NXF_HOME=${OUTPUT_DIR}
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
     -profile "sge_nih_trek" \
     --file_paths_10x "${STUDY_DIR}/file_paths_10x.tsv" \
     --file_metadata "${STUDY_DIR}/conditions_clinical_cr_metrics.tsv" \
     --file_sample_qc "${STUDY_DIR}/preprocessing_configs/params-sample_qc-doublets_parker_cutoffs.yml" \
     --genes_exclude_hvg "${STUDY_DIR}/../data-variable_gene_filter_ribo_mt.tsv" \
     --genes_score "${STUDY_DIR}/../data-gene_scores.tsv" \
     --output_dir "${OUTPUT_DIR}" \
     -params-file "${STUDY_DIR}/preprocessing_configs/params.yml" \
     --run_multiplet \
     -with-report \
     -with-trace \
     -with-timeline \
     -with-dag flowchart.png \
     -resume
