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


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object with principal components. Plot the PC loadings \
            principal component, plot top genes ranked by their contribution.
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
        '--num_pcs}',
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

    cell_type = np.unique(df['cell_type'])
    for i in cell_type:
        marker_genes = df.loc[df['cell_type'] == i]['hgnc_symbol']
        marker_genes_found = adata.var['gene_symbols'][
            adata.var['gene_symbols'].isin(marker_genes)
        ]
        # If there are loads of markers, just take the top 100
        if marker_genes_found.size > 100:
            marker_genes_found = marker_genes_found[0:100]

        # Clean up cell id names
        i_out = re.sub(r'\W+', '', i)  # Strip all non alphanumeric characters
        print(i, i_out)

        # Dotplots
        _ = sc.pl.dotplot(
            adata=adata,
            var_names=marker_genes_found,
            groupby='cluster',
            gene_symbols='gene_symbols',
            dendrogram=True,
            show=False,
            use_raw=False,
            log=False,
            color_map='Blues',
            save='-{}-{}-{}.png'.format(
                out_file_base,
                i_out,
                data_scale
            )
        )
        _ = sc.pl.dotplot(
            adata=adata,
            var_names=marker_genes_found,
            groupby='cluster',
            gene_symbols='gene_symbols',
            dendrogram=True,
            show=False,
            standard_scale='var',  # Scale color between 0 and 1
            use_raw=False,
            color_map='Blues',
            save='-{}-{}-{}_standardized.png'.format(
                out_file_base,
                i_out,
                data_scale
            )
        )


if __name__ == '__main__':
    main()
