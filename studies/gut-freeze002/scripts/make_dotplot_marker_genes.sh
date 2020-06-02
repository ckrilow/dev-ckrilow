#!/bin/sh

# this shell will kill itself if it takes over __ kb of virtual memory
# note: 8388608 kb = 8 gb
# ulimit -v 8388608


# --------------------------- Get and check user input -------------------------
dotplot_script=$1 # Example: dotplot-marker_database.py
adata_file=$2 # Example: "adata-normalized_pca-bbknn-umap-clustered.h5ad"
list_of_marker_files=${3-"../data-dotplot_marker_genes.tsv"}

# this will require file1 to be present
if [ -z ${dotplot_script} ]; then
    echo "ERROR: dotplot_script required"
    exit
fi
if [ -z ${adata_file} ]; then
    echo "ERROR: adata_file required"
    exit
fi
# ------------------------------------------------------------------------------

adata_file=$(realpath ${adata_file})
root=$(pwd)
for i in $(cat ${list_of_marker_files}); do
    cd ${root}
    # Make the submission command
    cmd="python ${dotplot_script} --h5_anndata ${adata_file}";
    cmd="${cmd} --markers_database ${i}";
    echo ${cmd};

    # Print the out dir
    out_dir=$(echo ${i} | awk -F/ '{print $(NF-1)}')
    echo ${out_dir}
    mkdir -p ${out_dir}
    cd ${out_dir}

    # Run the command
    bsub -q normal -G team152 -g /${USER} \
        -M 38192 -R "select[mem>38192] rusage[mem=38192]" \
        -o "out.out" -e "err.err" \
        "${cmd}"
done
