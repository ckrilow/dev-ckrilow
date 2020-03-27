
# Description

The methods used in this module are described in `docs/methods.pdf`.

Below is the structure of the results directory.
```bash
nf-qc_normalize_cluster
├── normalization_001::description_of_params
│   ├── [files: data]
│   ├── reduced_dims-pca::description_of_params
│   │   ├── [files: data]
│   │   ├── [plots: umap]
│   │   ├── cluster_001::description_of_params
│   │   │   ├── [files: data,cluster_marker_genes]
│   │   │   └── [plots: umap]
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
* `scanpy_merge-dev.py`: Implement per sample filtering for scanpy merge.
* `scanpy_normalize_pca.py`: Enable regress out.
* `scanpy_cluster.py`: Choose proper data normalization for marker identification and plots (right now I think it uses the same data used for dimensionality reduction).
* `scanpy_cluster.py`: Currently for clustering, we can change method (leiden or louvain), resolution, and n_pcs. Are there other paremeters that need to be scaled over?


# Enhancement list

* Check phenotypes against predicted sex from gene expression.
* Add basic QC plots - try to do this in R from anndata frame?
* Scrublet functionality + add to metadata + cluster distributions
* Gene scores + add to metadata
* Add marker gene AUC like here http://www.nxn.se/valent/2018/3/5/actionable-scrna-seq-clusters
* Add summary ARI and LISI metrics computed over a list of many different cluster annotations?
* Add tSNE plots - rapid plots with OpenTSNE?


# Quickstart

Quickstart for deploying this pipeline on a local cluster.


## 1. Set up the environment

Install most of the packages via conda:
```bash
# The repo directory.
REPO_MODULE="${HOME}/repo/path/to/this/module"

# Install environment using Conda.
conda env create --name sc_qc_normalize_cluster --file ${REPO_MODULE}/env/environment.yml

# Activate the new Conda environment.
source activate sc_qc_normalize_cluster

# To update environment file:
#conda env export --no-builds | grep -v prefix | grep -v name > environment.yml
```

Conda does not work well for R packages. Install R packages via renv:
```bash
TODO
```


## 2. Prepare the input files (TODO)

The pipeline takes as input:
1. File1.
2. File2.

Prepare the File1.
```bash
python "${REPO_MODULE}/scripts/prepare_file1_script.py" \
    --input_file input-file1.tsv > input-file1.tsv
```

Edit the config files for the nextflow anlaysis. An example can be found in `example_runtime_setup/`.


## 3. Set up and run Nextflow

Run Nexflow locally (NOTE: if running on a vm you may need to set `export QT_QPA_PLATFORM='offscreen'` for scanpy as described [here](https://github.com/ipython/ipython/issues/10627)):
```bash
# NOTE: all input file paths should be full paths.
nextflow run "${REPO_MODULE}/main.nf" \
    -profile "local" \
    --file_paths_10x "/home/ubuntu/studies/TaylorDL/freeze001-dev/final_samples2.tsv" \
    --file_metadata "/home/ubuntu/repo/scrna_cellranger/sync_status/samples_metainfo.tsv" \
    -resume
```


Run Nextflow using LSF on a compute cluster. More on bgroups [here](https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsb.params.default_jobgroup.5.html).:
```bash
# Boot up tmux session.
tmux new -s nf

# Log into an interactive session.
bgadd "/${USER}/logins"
bsub -q normal -G team152 -g /lt9/logins -Is -XF -M 8192 -R "select[mem>8192] rusage[mem=8192]" /bin/bash

# Activate the Conda environment (inherited by subsequent jobs).
source activate sc_qc_normalize_cluster

# Set up a group to submit jobs to (export a default -g parameter).
bgadd -L 500 "/${USER}/nf"
export LSB_DEFAULT_JOBGROUP="/${USER}/nf"
# Depending on LSF setup, you may want to export a default -G parameter.
export LSB_DEFAULTGROUP="team152"
# NOTE: by setting the above flags, all of the nextflow lsf jobs will have these flags set.

# Change to a temporary runtime dir on the node:
mkdir -p "/tmp/${USER}/nf"
cd "/tmp/${USER}/nf"

# NOTE: if you want to resume a previous workflow, add -resume to the flag
nextflow run "${REPO_MODULE}/main.nf" \
    -profile "lsf" \
    --file_paths_10x "/home/ubuntu/studies/TaylorDL/freeze001-dev/final_samples2.tsv" \
    --file_metadata "/home/ubuntu/repo/scrna_cellranger/sync_status/samples_metainfo.tsv"
```


Example of how one may sync results:
```bash
NF_OUT_DIR="/path/to/out_dir"
rsync -am --include='*.png' --include='*/' --exclude='*' my_cluster_ssh:${NF_OUT_DIR} .
rsync -am --include='*.png' --include='*/' --exclude='*' my_cluster_ssh:${NF_OUT_DIR} .
```
