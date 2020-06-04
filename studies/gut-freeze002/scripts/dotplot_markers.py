#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-05-29'
__version__ = '0.0.1'

import argparse
import os
import scanpy as sc
import pandas as pd


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and markers from database. Plot expression \
            expression of markers in each cluster.
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
        '--markers_csv',
        action='store',
        dest='markers_csv',
        default='none',
        help='Markers to plot. Must have the following columns: cell_type, \
            hgnc_symbol.'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    # sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    sc.set_figure_params(dpi_save=300)

    # Get the out file base.
    out_file_base = os.path.basename(options.markers_csv).rstrip('csv').rstrip('.')

    # Read in the data
    adata = sc.read_h5ad(filename=options.h5)

    # Read in the marker database file
    df = pd.read_table(options.markers_csv)
    df = df['hgnc_symbol']

    marker_genes_found = adata.var['gene_symbols'][
            adata.var['gene_symbols'].isin(df)
        ]

    # Dotplots
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
            save='{}_dotplot.png'.format(out_file_base)
    )


if __name__ == '__main__':
    main()
