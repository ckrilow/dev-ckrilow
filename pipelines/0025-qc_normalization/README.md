
# Description

The methods used in this module are described in `docs/methods.pdf`.


# TODO list

* Add `docs/methods.pdf` file.
* Add brief description of module.
* `scanpy_merge-dev.py`: Implement per sample filtering for scanpy merge.
* `scanpy_normalize_pca.py`: Fix highly variable gene gene detection in pca.
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
# the repo directory
REPO_MODULE="${HOME}/repo/path/to/this/module"

# install environment using Conda
conda env create --name sc_qc_normalization_cluster --file ${REPO_MODULE}/environment.yml

# activate the new Conda environment
source activate sc_qc_normalization_cluster

# to update environment file:
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


## 3. Set up and run Nextflow

Change to a temporary runtime dir:
```bash
mkdir -p "/tmp/${USER}/nf"
cd "/tmp/${USER}/nf"
```


Run Nexflow locally:
```bash
# NOTE: all input file paths should be full paths.
nextflow run "${REPO_MODULE}/main.nf" \
    -profile "local" \
    --file_paths_10x "/home/ubuntu/studies/TaylorDL/freeze001-dev/final_samples2.tsv" \
    --file_metadata "/home/ubuntu/repo/scrna_cellranger/sync_status/samples_metainfo.tsv" \
    -resume
```


Run Nextflow using LSF on a local cluster. More on bgroups [here](https://www.ibm.com/support/knowledgecenter/SSETD4_9.1.3/lsf_config_ref/lsb.params.default_jobgroup.5.html):
```bash
# Set up a group to submit jobs to (export a default -g parameter).
bgadd -L 500"/${USER}/nf"
export LSB_DEFAULT_JOBGROUP="/${USER}/nf"

# Depending on LSF setup, you may want to export a default -G parameter.
#export LSB_DEFAULTGROUP="team152"

# NOTE: if you want to resume a previous workflow, add -resume to the flags
nextflow run "${REPO_MODULE}/main.nf" \
    -profile "local" \
    --file_paths_10x "/home/ubuntu/studies/TaylorDL/freeze001-dev/final_samples2.tsv" \
    --file_metadata "/home/ubuntu/repo/scrna_cellranger/sync_status/samples_metainfo.tsv"
```


Example of how one may sync results:
```bash
NF_OUT_DIR="/path/to/out_dir"
rsync -am --include='*.png' --include='*/' --exclude='*' my_cluster:${NF_OUT_DIR} .
rsync -am --include='*.png' --include='*/' --exclude='*' my_cluster:${NF_OUT_DIR} .
```
