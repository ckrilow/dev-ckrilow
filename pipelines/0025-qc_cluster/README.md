
# Description

The methods used in this module are described in `docs/methods.pdf`. TODO: `docs/methods.pdf`

Below is the structure of the results directory. The values that will be listed in `description_of_params` within the directory structure correspond to the various parameters one can set. An example of a paramters file is found in `example_runtime_setup/params.yml`.
```bash
nf-qc_cluster
├── normalization_001::description_of_params
│   ├── [files: data]
│   ├── reduced_dims-pca::description_of_params
│   │   ├── [files: data]
│   │   ├── [plots: umap]
│   │   ├── cluster_001::description_of_params
│   │   │   ├── [files: data,clusters]
│   │   │   ├── [plots: umap]
│   │   │   ├── cluster_markers_001::description_of_params
│   │   │   │   ├── [files: cluster_marker_genes]
│   │   │   │   └── [plots: marker_genes,marker_genes_dotplot]
│   │   │   ├── cluster_markers_002::description_of_params
│   │   │   ... etc. ...
│   │   ├── cluster_002::description_of_params
│   │   ... etc. ...
│   ├── reduced_dims-harmony_001::description_of_params
│   ├── reduced_dims-harmony_002::description_of_params
│   ... etc. ...
├── normalization_002::description_of_norm_params
... etc. ...
└── adata.h5  # concatenated single cell data with no normalization
```


# TODO list

* Add `docs/methods.pdf` file.
* Add brief description of module.


# Enhancement list

* `scanpy_merge-dev.py`: If it were important to have a per sample filter, merge could be re-designed to accommodate this.
* `scanpy_cluster.py`: Currently for clustering, we can change method (leiden or louvain), resolution, and n_pcs. Are there other parameters that need to be scaled over?
* Check phenotypes against predicted sex from gene expression.
* Add basic QC plots - try to do this in R from anndata frame?
* Scrublet functionality + add to metadata + cluster distributions
* Gene scores + add to metadata
* Add marker gene AUC like here http://www.nxn.se/valent/2018/3/5/actionable-scrna-seq-clusters
* Add summary ARI and LISI metrics computed over a list of many different cluster annotations?
* Add tSNE plots - rapid plots with OpenTSNE?
* Calculate marker genes with diffxpy or logreg?


# Quickstart

Quickstart for deploying this pipeline locally and on a high performance compute cluster.


## 1. Set up the environment

Install the required packages via conda:
```bash
# The repo directory.
REPO_MODULE="${HOME}/repo/path/to/this/pipeline"

# Install environment using Conda.
conda env create --name sc_qc_cluster --file ${REPO_MODULE}/env/environment.yml

# Activate the new Conda environment.
source activate sc_qc_cluster

# To update environment file:
#conda env export --no-builds | grep -v prefix | grep -v name > environment.yml
```


## 2. Prepare the input files

Generate and/or edit input files for the pipeline.

The pipeline takes as input:
1. **file_paths_10x**:  Tab-delimited file containing experiment_id and data_path_10x_format columns. Reqired.
2. **file_metadata**:  Tab-delimited file containing sample metadata. Reqired.
3. **file_sample_qc**:  YAML file containing sample qc and filtering parameters. Required. NOTE: in the example config file, this is part of the YAML file for item 4.
4. **params-file**:  YAML file containing analysis parameters. Optional.

Examples of all of these files can be found in `example_runtime_setup/`.


## 3. Set up and run Nextflow

Run Nexflow locally (NOTE: if running on a virtual machine you may need to set `export QT_QPA_PLATFORM="offscreen"` for scanpy as described [here](https://github.com/ipython/ipython/issues/10627)):
```bash
# Boot up tmux session.
tmux new -s nf

# NOTE: All input file paths should be full paths.
nextflow run "${REPO_MODULE}/main.nf" \
    -profile "local" \
    --file_paths_10x "${REPO_MODULE}/example_runtime_setup/file_paths_10x.tsv" \
    --file_metadata "${REPO_MODULE}/example_runtime_setup/file_metadata.tsv" \
    -params-file "${REPO_MODULE}/example_runtime_setup/params.yml"
```


Run Nextflow using LSF on a compute cluster. More on bgroups [here](https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsb.params.default_jobgroup.5.html).:
```bash
# Set the results directory.
RESULTS_DIR="/path/to/results/dir"
mkdir -p "${RESULTS_DIR}"

# Boot up tmux session.
tmux new -s nf

# Log into an interactive session.
# NOTE: Here we set the -G parameter due to our institute's LSF configuration.
bgadd "/${USER}/logins"
bsub -q normal -G team152 -g /${USER}/logins -Is -XF -M 8192 -R "select[mem>8192] rusage[mem=8192]" /bin/bash
# NOTE: If you are running over many cells, you may need to start an
# interactive session on a queue that allows long jobs
#bsub -q long -G team152 -g /${USER}/logins -Is -XF -M 18192 -R "select[mem>18192] rusage[mem=18192]" /bin/bash

# Activate the Conda environment (inherited by subsequent jobs).
conda activate sc_qc_cluster

# Set up a group to submit jobs to (export a default -g parameter).
bgadd -L 500 "/${USER}/nf"
export LSB_DEFAULT_JOBGROUP="/${USER}/nf"
# Depending on LSF setup, you may want to export a default -G parameter.
export LSB_DEFAULTGROUP="team152"
# NOTE: By setting the above flags, all of the nextflow LSF jobs will have
# these flags set.

# Settings for scanpy (see note above).
export QT_QPA_PLATFORM="offscreen"

# Change to a temporary runtime dir on the node. In this demo, we will change
# to the same execution directory.
cd ${RESULTS_DIR}

# Remove old logs and nextflow output (if one previously ran nextflow in this
# dir).
rm -r *html;
rm .nextflow.log*;

# NOTE: If you want to resume a previous workflow, add -resume to the flag.
nextflow run "${REPO_MODULE}/main.nf" \
    -profile "lsf" \
    --file_paths_10x "${REPO_MODULE}/example_runtime_setup/file_paths_10x.tsv" \
    --file_metadata "${REPO_MODULE}/example_runtime_setup/file_metadata.tsv" \
    --file_sample_qc "${REPO_MODULE}/example_runtime_setup/params.yml" \
    --output_dir "${RESULTS_DIR}" \
    -params-file "${REPO_MODULE}/example_runtime_setup/params.yml" \
    -with-report "nf_report.html" \
    -resume

# NOTE: If you would like to see the ongoing processes, look at the log files.
cat .nextflow.log
```


Example of how one may sync results:
```bash
NF_OUT_DIR="/path/to/out_dir"
rsync -am --include="*.png" --include="*/" --exclude="*" my_cluster_ssh:${NF_OUT_DIR} .
rsync -am --include="*.png" --include="*/" --exclude="*" my_cluster_ssh:${NF_OUT_DIR} .
```
