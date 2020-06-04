
# Description

Basic usage of dotplot_markers.py

# How to run

Prepare .csv file with markers of interest (see "example_markers.csv"). Note
only one column is requried "hgnc_symbol". Then run the .py file:


python ./dotplot_markers.py --h5_anndata "adata-normalized_pca-bbknn-umap-clustered.h5ad" --markers_csv "example_markers.csv"

Assumption is that all three files (.py, .h5ad and .csv are in the same location).

# Requirements

To run the script you need to have installed following packages:
* argparse
* os
* scanpy
* pandas

Install them using i.e. pip install argparse.
