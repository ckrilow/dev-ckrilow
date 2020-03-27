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
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-pc', '--tsv_pcs',
        action='store',
        dest='pc',
        default='',
        help='Tab-delimited file of PCs for each cell. First column is\
            cell_barcode. Subsequent columns are PCs. If "", uses pca\
            slot in AnnData.\
            (default: "")'
    )

    parser.add_argument(
        '-cm', '--cluster_method',
        action='store',
        dest='cm',
        default='leiden',
        help='Clustering method. Valid options: [leiden|louvain].\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-npc', '--number_pcs',
        action='store',
        dest='npc',
        default=0,
        type=int,
        help='Number of PCs to use.\
            (default: maximum number in tsv_pcs file)'
    )

    parser.add_argument(
        '-r', '--resolution',
        action='store',
        dest='r',
        default=1.0,
        type=float,
        help='Resolution.\
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
    verbose = True

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-{}-clustered'.format(
            os.path.basename(options.h5.rstrip('.h5')),
            os.path.basename(options.pc.rstrip('.tsv.gz'))
        )

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)

    # Load the PCs.
    if options.pc == '':
        df_pca = pd.DataFrame(
            data=adata.obsm['X_pca'],
            index=adata.obs.index,
            columns=[
                'PC{}'.format(x) for x in
                range(1, adata.obsm['X_pca'].shape[1]+1)
            ]
        )
    else:
        df_pca = pd.read_csv(options.pc, sep='\t', index_col='cell_barcode')

    # Check that nPCs is valid.
    n_pcs = options.npc
    if n_pcs == 0:
        n_pcs = len(df_pca.columns)
    elif n_pcs > len(df_pca.columns):
        raise Exception(
            '--number_pcs ({}) is > than n_pcs in --tsv_pcs ({}).'.format(
                n_pcs,
                len(df_pca.columns)
            )
        )
    if verbose:
        print('Using {} PCs.'.format(n_pcs))

    # Add the reduced dimensions to the AnnData object.
    adata.obsm['X_pca'] = df_pca.loc[adata.obs.index, :].values.copy()

    # Calculate neighbors for on the specified PCs.
    sc.pp.neighbors(
        adata,
        use_rep='X_pca',
        # n_neighbors=10,
        n_pcs=n_pcs,
        copy=False
    )

    # Run the clustering, choosing either leiden or louvain algorithms
    cluster_method = options.cm
    cluster_resolution = options.r
    if cluster_method == 'leiden':
        sc.tl.leiden(
            adata,
            resolution=cluster_resolution,
            key_added=cluster_method,
            copy=False
        )
    elif cluster_method == 'louvain':
        sc.tl.louvain(
            adata,
            flavor='vtraag',
            resolution=cluster_resolution,
            key_added=cluster_method,
            copy=False
        )
    else:
        raise Exception(
            'Invalid --cluster_method: {}.'.format(
                cluster_method
            )
        )
    # Also save the clusters to the same spot so we know where they will be.
    adata.uns['cluster'] = adata.uns[cluster_method]
    adata.uns['cluster']['params']['method'] = cluster_method
    adata.obs['cluster'] = adata.obs[cluster_method]

    # Print the final number of clustered discrovered
    if verbose:
        print('{} clusters identified'.format(
            adata.obs[cluster_method].drop_duplicates().shape[0]
        ))

    # Save the clustered data to a data frame.
    cell_clustering_df = adata.obs[[cluster_method]].copy()
    cell_clustering_df.columns = ['cluster']
    cell_clustering_df['cluster_method'] = cluster_method
    cell_clustering_df['cluster_resolution'] = cluster_resolution
    cell_clustering_df.to_csv(
        '{}.tsv.gz'.format(out_file_base),
        sep='\t',
        index=True,
        quoting=csv.QUOTE_NONNUMERIC,
        index_label='cell_barcode',
        na_rep='',
        compression='gzip'
    )

    adata.write('{}.h5'.format(out_file_base), compression='gzip')

    if verbose:
        print('Finished clustering and saved clustered AnnData.')

    # TODO: Check what data is being used here: raw, lognorm, or norm+scaled?
    #
    # Identify cell type makers.
    # Wilcoxon is recommended over t-test in this paper:
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

    # Other options for cell type marker identification:
    # MAST, limma, DESeq2, diffxpy, logreg

    # Rank genes by logistic regression (multi-variate appraoch), suggested by:
    # https://doi.org/10.1101/258566
    # NOTE: adata.uns['rank_genes_groups'].keys() = ['params',
    #       'scores', 'names']
    # sc.tl.rank_genes_groups(
    #     adata,
    #     groupby=cluster_method,
    #     groups='all',
    #     reference='rest',
    #     method='logreg'
    # )

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
