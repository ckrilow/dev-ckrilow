#!/usr/bin/env python


__author__ = 'Tobi Alegbe'
__date__ = '2020-06-30'
__version__ = '0.0.1'

import argparse
import os
import re
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec

# Nice large palette.
COLORS_LARGE_PALLETE = [
    '#0F4A9C', '#3F84AA', '#C9EBFB', '#8DB5CE', '#C594BF', '#DFCDE4',
    '#B51D8D', '#6f347a', '#683612', '#B3793B', '#357A6F', '#989898',
    '#CE778D', '#7F6874', '#E09D37', '#FACB12', '#2B6823', '#A0CC47',
    '#77783C', '#EF4E22', '#AF1F26'
]

def save_pc_fig(
     adata,
    dict__umap_dim_and_params,
    out_file_base,
    color_var,
    colors_quantitative=True,
    colors_large_palette=COLORS_LARGE_PALLETE,
    drop_legend=-1
):
    """
    Create and save a figure containing the cells plotted in the space of \
    two principal components.
    """

def save_pc_genes_fig(
     adata,
    dict__umap_dim_and_params,
    out_file_base,
    color_var,
    colors_quantitative=True,
    colors_large_palette=COLORS_LARGE_PALLETE,
    drop_legend=-1
):
    """
    Create and save a figure containing the top genes that contribute to a \
    principal component.
    """
    
    # Plot the top genes contributing to each principal component
    adata_temp = adata.copy()
    adata_temp.var_names = adata_temp.var.gene_symbols
    pcs_to_plot = ','.join(map(str,list(range(1,num_PCs))))

    sc.pl.pca_loadings(
        adata=adata_temp,
        include_lowest=False,
        components = pcs_to_plot
        )

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object with principal components. Plot the cells in \
            principal component coordinates and plot top genes in each component\
            ranked by their contribution.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '--num_pcs',
        action='store',
        dest='num_PCs',
        default='20',
        help='Number of principal components to plot'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: None)'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    # sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    sc.set_figure_params(dpi_save=300)

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.'))
        )

    # Read in the data
    adata = sc.read_h5ad(filename=options.h5)

    # Plot the cells in the principal component space
    pc_pairs_to_plot = ["{},{}".format(i,i+1) for i in range(1,options.num_PCs,2)]

    sc.pl.pca(
        adata=adata,
        color=color_var,
        palette=color_palette,
        components=pc_pairs_to_plot,
        alpha=0.4
        )


    fig.savefig(
        '{}-{}.png'.format(out_file_base, color_var),
        dpi=300,
        bbox_inches='tight'
    )


if __name__ == '__main__':
    main()
