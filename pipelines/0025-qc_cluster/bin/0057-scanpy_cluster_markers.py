#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scanpy as sc
import csv


def nd(arr):
    """Make nd array."""
    return np.asarray(arr).reshape(-1)


def dotplot(adata, markers):
    """Customize dotplot function."""
    adata_raw = adata.raw.to_adata()
    cluter_ids = np.unique(
        adata_raw.obs.cluster.values.astype(int)
    ).astype(str)

    features = adata_raw.var.gene_symbols.values
    midx = [np.where(i == features)[0][0] for i in markers]

    # for each cluster for each gene get two things
    # 1 percent of cells in the cluster expressing that gene
    # 2 average expression of that gene (for cells that are expressing it)
    per = np.zeros((len(cluter_ids), len(markers)))
    avg = np.zeros((len(cluter_ids), len(markers)))
    mtx = adata_raw.X  # .todense()

    for cn, c in enumerate(cluter_ids):
        tmp_mtx = mtx[adata_raw.obs.cluster.values == c]
        sub_mtx = tmp_mtx[:, midx]
        avg[cn] = nd(sub_mtx.mean(axis=0))
        per[cn] = (sub_mtx > 0).sum(axis=0) / sub_mtx.shape[0]

    fig, ax = plt.subplots(
        figsize=(len(markers) * 0.35, len(cluter_ids) * 0.3 + 1)
    )
    xidx = np.arange(len(markers))
    yidx = np.arange(len(cluter_ids))

    xlabels = markers
    ylabels = cluter_ids

    X, Y = np.meshgrid(xidx, yidx)

    for dn, d in enumerate(per):
        # a = ax.scatter(X[dn], Y[dn], s=d*500+10, c=avg[dn], cmap='Blues')
        a = ax.scatter(X[dn], Y[dn], s=d*100, c=avg[dn], cmap='Blues')

    ax.set_xticks(xidx)
    ax.set_yticks(yidx)

    ax.set_xticklabels(xlabels, rotation=90, ha='center')
    ax.set_yticklabels(ylabels)

    ax.set_xlabel('Gene')
    ax.set_ylabel('Cluster')
    ax.figure.colorbar(a, ax=ax, label='$ln(CPM+1)$')

    handles = [
        mpl.lines.Line2D([
            0], [0], marker='o', color='w', label='  0%',
            markerfacecolor='black', markersize=7
        ),
        mpl.lines.Line2D(
            [0], [0], marker='o', color='w', label=' 25%',
            markerfacecolor='black', markersize=10
        ),
        mpl.lines.Line2D(
            [0], [0], marker='o', color='w', label=' 50%',
            markerfacecolor='black', markersize=12
        ),
        mpl.lines.Line2D(
            [0], [0], marker='o', color='w', label=' 75%',
            markerfacecolor='black', markersize=13.5
        ),
        mpl.lines.Line2D(
            [0], [0], marker='o', color='w', label='100%',
            markerfacecolor='black', markersize=17
        )
    ]
    ax.legend(
        handles=handles,
        loc='center left',
        bbox_to_anchor=(1.15, 0.5),
        title='Percent of cells'
    )

    # plt.show()

    return(fig, ax)


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

    print(
        'WARNING: All functions in this script set use_raw=True,',
        'assuming that adata.raw.to_adata stores ln(CPM+1) normalized data.'
    )

    # NOTE: You should be using the ln(CPM+1) data here. By default these
    #       functions use the .raw attribute of AnnData if present which is
    #       assumed to be ln(CPM+1).

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
            use_raw=True,
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
            use_raw=True,
            method='logreg',
            n_genes=100,
            max_iter=5000  # passed to sklearn.linear_model.LogisticRegression
        )
    else:
        # Other options for cell type marker identification in scanpy:
        # MAST, limma, DESeq2, diffxpy, logreg
        raise Exception('Method not implemented')

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
    # NOTE: Because rank_genes_groups was run on the ln(CPM+1) data,
    #       we must be sure to use the gene symbols from that data since the
    #       data in adata.X may contain fewer genes, for instance if the
    #       matrix was filtered before scaling.
    marker_df = marker_df.set_index('ensembl_gene_id', inplace=False)
    marker_df = marker_df.join(
        adata.raw.var[['gene_symbols']],
        how='left'
    )
    if (np.invert(marker_df.gene_symbols.notnull()).sum() > 0):
        filt = np.invert(marker_df.gene_symbols.notnull())
        print(marker_df.loc[filt, :])
        raise Exception('Missing gene_symbols in marker_df.')
    elif (np.invert(marker_df.gene_symbols.notna()).sum() > 0):
        filt = np.invert(marker_df.gene_symbols.notna())
        print(marker_df.loc[filt, :])
        raise Exception('Missing gene_symbols in marker_df.')

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

    # Plot cell type makers.
    # Annoyingly, prefix hardcoded as rank_genes_groups_<cluster_id>.
    _ = sc.pl.rank_genes_groups(
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
    if (np.invert(marker_df_plt.gene_symbols.notnull()).sum() > 0):
        filt = np.invert(marker_df_plt.gene_symbols.notna())
        print(marker_df_plt.loc[filt, :])
        raise Exception('Missing gene_symbols in marker_df_plt.')
    elif (np.invert(marker_df_plt.gene_symbols.notna()).sum() > 0):
        filt = np.invert(marker_df_plt.gene_symbols.notna())
        print(marker_df_plt.loc[filt, :])
        raise Exception('Missing gene_symbols in marker_df_plt.')

    # NOTE: You should be using the ln(CPM+1) data here. "$ln(CPM+1)$"
    # Plot cell type markers in dotplot.
    fig, ax = dotplot(
        adata,
        markers=marker_df_plt['gene_symbols'].values
    )
    fig.savefig(
        'dotplot-{}.pdf'.format(out_file_base),
        dpi=300,
        bbox_inches='tight'
    )
    # Annoyingly, the prefix is hardcoded as dotplot. Should save to axes
    _ = sc.pl.dotplot(
        adata,
        marker_df_plt['ensembl_gene_id'].to_list(),
        groupby='cluster',
        dendrogram=True,
        # gene_symbols='gene_symbols',  # causes error for some reason
        use_raw=True,
        show=False,
        color_map='Blues',
        save='_ensembl-{}.pdf'.format(out_file_base)
    )
    _ = sc.pl.dotplot(
        adata,
        marker_df_plt['ensembl_gene_id'].to_list(),
        groupby='cluster',
        dendrogram=True,
        standard_scale='var',  # Scale color between 0 and 1
        use_raw=True,
        show=False,
        color_map='Blues',
        save='_ensembl-{}-standardized.pdf'.format(out_file_base)
    )
    _ = sc.pl.heatmap(
        adata,
        marker_df_plt['ensembl_gene_id'].to_list(),
        groupby='cluster',
        dendrogram=True,
        use_raw=True,
        show=False,
        save='_ensembl-{}-standardized.pdf'.format(out_file_base)
    )
    _ = sc.pl.heatmap(
        adata,
        marker_df_plt['ensembl_gene_id'].to_list(),
        groupby='cluster',
        dendrogram=True,
        standard_scale='var',  # Scale color between 0 and 1
        use_raw=True,
        show=False,
        save='_ensembl-{}-standardized.pdf'.format(out_file_base)
    )


if __name__ == '__main__':
    main()
