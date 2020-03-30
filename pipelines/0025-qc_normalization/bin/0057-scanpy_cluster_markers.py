#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import pandas as pd
import scanpy as sc
import csv


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and PCs file. Clusters the data. Saves an
            AnnData object with clusters in the clusters slot, a clusters
            file, and QC plots.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file where clusters have been saved to cluster slot.'
    )

    parser.add_argument(
        '-rgm', '--rank_genes_method',
        action='store',
        dest='rgm',
        default='wilcoxon',
        help='Method used to rank marker genes. Valid options:\
            [wilcoxon|logreg].\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=4,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <h5_anndata>-<tsv_pcs>-clustered)'
    )

    options = parser.parse_args()

    # Fixed settings.
    # verbose = True

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}'.format(
            os.path.basename(options.h5.rstrip('.h5'))
        )

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)

    # TODO: Check what data is being used here: raw, lognorm, or norm+scaled?
    #
    # Identify cell type makers.
    if options.rgm == 'wilcoxon':
        # If considering simple Wilcoxon vs t-test, choose Wilcoxon as
        # recommended in this paper:
        # https://www.nature.com/articles/nmeth.4612.
        # NOTE: adata.uns['rank_genes_groups'].keys() = ['params',
        #       'scores', 'names', 'logfoldchanges', 'pvals', 'pvals_adj']
        sc.tl.rank_genes_groups(
            adata,
            groupby='cluster',
            groups='all',
            reference='rest',
            method='wilcoxon',
            n_genes=100,
            corr_method='bonferroni'
        )
    elif options.rgm == 'logreg':
        # Rank genes by logistic regression (multi-variate appraoch),
        # suggested by: https://doi.org/10.1101/258566
        # NOTE: adata.uns['rank_genes_groups'].keys() = ['params',
        #       'scores', 'names']
        sc.tl.rank_genes_groups(
            adata,
            groupby='cluster',
            groups='all',
            reference='rest',
            method='logreg',
            n_genes=100
        )
    else:
        # Other options for cell type marker identification in scanpy:
        # MAST, limma, DESeq2, diffxpy, logreg
        raise Exception(
            'Method not implemented'
        )

    # Save the ranks.
    results_dict = dict()
    for cluster_i in adata.uns['rank_genes_groups']['names'].dtype.names:
        # print(cluster_i)
        # Get keys that we want from the dataframe.
        data_keys = list(
            set(['names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']) &
            set(adata.uns['rank_genes_groups'].keys())
        )
        # Build a table using these keys.
        key_i = data_keys.pop()
        results_dict[cluster_i] = pd.DataFrame(
            row[cluster_i] for row in adata.uns['rank_genes_groups'][key_i]
        )
        results_dict[cluster_i].columns = [key_i]
        for key_i in data_keys:
            results_dict[cluster_i][key_i] = [
                row[cluster_i] for row in adata.uns['rank_genes_groups'][key_i]
            ]
        results_dict[cluster_i]['cluster'] = cluster_i
    marker_df = pd.concat(results_dict, ignore_index=True)

    # Clean up naming.
    marker_df = marker_df.rename(
        columns={'names': 'ensembl_gene_id'},
        inplace=False
    )

    # Add gene_symbols.
    marker_df = marker_df.set_index('ensembl_gene_id', inplace=False)
    marker_df = marker_df.join(
        adata.var[['gene_symbols']],
        how='left'
    )
    marker_df = marker_df.reset_index(drop=False)
    marker_df = marker_df.rename(
        columns={'index': 'ensembl_gene_id'},
        inplace=False
    )

    # Save the marker dataframe.
    marker_df.to_csv(
        '{}-cluster_markers.tsv.gz'.format(out_file_base),
        sep='\t',
        index=True,
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression='gzip'
    )

    # TODO: Check what data is being used here: raw, lognorm, or norm+scaled?
    # Plot cell type makers.
    # Annoyingly, prefix hardcoded as rank_genes_groups_<cluster_id>.
    sc.pl.rank_genes_groups(
        adata,
        gene_symbols='gene_symbols',
        n_genes=25,
        sharey=False,
        show=False,
        save='-{}.pdf'.format(out_file_base)
    )

    # Plot cell type markers in dotplot.
    # Annoyingly, prefix hardcoded as dotplot.
    # sc.pl.rank_genes_groups_dotplot(
    #     adata,
    #     n_genes=25,
    #     sharey=False,
    #     show=False,
    #     save='-{}.pdf'.format(out_file_base)
    # )

    # Sort by scores: same order as p-values except most methods return scores.
    marker_df = marker_df.sort_values(by=['scores'], ascending=False)
    # Make dataframe of the top 5 markers per cluster
    marker_df_plt = marker_df.groupby('cluster').head(3)

    # TODO: Check what data is being used here: raw, lognorm, or norm+scaled?
    # Plot cell type markers in dotplot.
    # Annoyingly, prefix hardcoded as dotplot.
    sc.pl.dotplot(
        adata,
        marker_df_plt['ensembl_gene_id'].to_list(),
        groupby='cluster',
        dendrogram=True,
        # gene_symbols='gene_symbols',
        use_raw=True,
        show=False,
        save='_ensembl-{}.pdf'.format(out_file_base)
    )


if __name__ == '__main__':
    main()
