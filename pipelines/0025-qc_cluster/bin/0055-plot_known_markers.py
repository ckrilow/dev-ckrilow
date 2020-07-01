#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-05-29'
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
        '--markers_database_tsv',
        action='store',
        dest='markers_database_tsv',
        default='none',
        help='Markers to plot. Must have the following columns: cell_type, \
            hgnc_symbol.'
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

    data_scale = 'log1p_cp10k'
    if data_scale == 'cp10k':  # Set X to cp10k
        adata.X = np.expm1(adata.layers['log1p_cp10k'])
    elif data_scale == 'log1p_cp10k':  # Set X to log1p_cp10k
        adata.X = adata.layers['log1p_cp10k']
    elif data_scale == 'counts':
        adata.X = adata.layers['counts']
        #adata = adata.raw.to_adata() # This is not counts.

    # Read in the marker database file
    df = pd.read_table(options.markers_database_tsv)

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
        # If there are loads of markers, just take the top 100
        if marker_genes_found.size > 100:
            marker_genes_found = marker_genes_found[0:100]

        # Clean up cell id names
        i_out = re.sub(r'\W+', '', i)  # Strip all non alphanumeric characters
        print(i, i_out)

        # Dotplots
        # Generate dendrogram using the marker genes... this will be used in the
        # below dotplots.
        # NOTE: With latest version of pandas, sc.tl.dendrogram throws an error.
        run_dendrogram = False
        if run_dendrogram:
            sc.tl.dendrogram(
                adata,
                groupby='cluster',
                use_rep='X_pca',
                var_names=marker_dict_plt,
                use_raw=False,
                cor_method='pearson',
                linkage_method='complete',
                optimal_ordering=True,
                inplace=True
            )
        
        _ = sc.pl.dotplot(
            adata=adata,
            var_names=marker_genes_found,
            groupby='cluster',
            gene_symbols='gene_symbols',
            dendrogram=run_dendrogram,
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
            dendrogram=run_dendrogram,
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

        # Plot heatmaps
        # _ = sc.pl.heatmap(
        #     adata=adata_raw,
        #     var_names=marker_genes_found,
        #     groupby='cluster',
        #     dendrogram=True,
        #     gene_symbols='gene_symbols',
        #     use_raw=False,
        #     show=False,
        #     save='-{}-{}.png'.format(out_file_base, i_out)
        # )
        # _ = sc.pl.heatmap(
        #     adata=adata_raw,
        #     var_names=marker_genes_found,
        #     groupby='cluster',
        #     dendrogram=True,
        #     gene_symbols='gene_symbols',
        #     standard_scale='var',  # Scale color between 0 and 1
        #     use_raw=False,
        #     show=False,
        #     save='-{}-{}-standardized.png'.format(out_file_base, i_out)
        # )


if __name__ == '__main__':
    main()
