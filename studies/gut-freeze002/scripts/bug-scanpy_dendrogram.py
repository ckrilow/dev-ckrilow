#!/usr/bin/env python

__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import os
import random
import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
matplotlib.use('Agg')

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


def main():

    # Use pbmc3k dataset
    adata = sc.datasets.pbmc3k()
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.log1p(adata)
    sc.pp.normalize_total(adata)
    sc.pp.highly_variable_genes(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    sc.tl.rank_genes_groups(adata, groupby='leiden')

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

    # Sort by scores: same order as p-values except most methods return scores.
    marker_df = marker_df.sort_values(by=['scores'], ascending=False)
    # Make dataframe of the top 3 markers per cluster
    marker_df_plt = marker_df.groupby('cluster').head(3)

    _ = sc.pl.dotplot(
        adata,
        var_names=marker_df_plt['names'],
        groupby='leiden',
        dendrogram=True,
        use_raw=False,
        show=False,
        color_map='Blues'
        #save='{}.png'.format('test')
    )


if __name__ == '__main__':
    main()
