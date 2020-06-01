#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-05-29'
__version__ = '0.0.1'

import argparse
import os
import scanpy as sc
import pandas as pd
import numpy as np


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
        '--markers_database',
        action='store',
        dest='markers_database',
        default='none',
        help='Markers from database to be dotplotted.'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    # sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    sc.set_figure_params(dpi_save=300)

    adata = sc.read_h5ad(filename=options.h5)

    df = pd.read_table(options.markers_database)

    if 'p_value_adj' in df.columns:
        df = df[['cell_type', 'hgnc_symbol', 'p_value_adj']]
        df = df.sort_values(
            by=['cell_type', 'p_value_adj'],
            ascending=[True, True]
        )
    else:
        df = df[['cell_type', 'hgnc_symbol']]

    cell_type = np.unique(df['cell_type'])
    for i in cell_type:
        marker_genes = df.loc[df['cell_type'] == i]['hgnc_symbol']
        marker_genes_found = adata.var['gene_symbols'][
            adata.var['gene_symbols'].isin(marker_genes)
        ]
        if marker_genes_found.size > 100:
            marker_genes_found = marker_genes_found[0:100]
        i_out = i.replace(' ', '_').replace('.', '_').replace('/', '-')
        print(i, i_out)
        sc.pl.dotplot(
            adata,
            marker_genes_found,
            groupby='cluster',
            gene_symbols='gene_symbols',
            dendrogram=True,
            save='-{}.png'.format(i_out)
        )


if __name__ == '__main__':
    main()
