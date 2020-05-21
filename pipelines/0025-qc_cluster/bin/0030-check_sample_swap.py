#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-05-20'
__version__ = '0.0.1'

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import plotnine as plt9
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (10, 10)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read anndata object. Read chrX and chrY genes, plot scatterplot \
            of mean expression of those signature across samples.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-Y', '--chrY_genes',
        required=True,
        help='CSV file.'
    )

    parser.add_argument(
        '-X', '--chrX_genes',
        required=True,
        help='CSV file.'
    )

    parser.add_argument(
        '-o', '--output.file',
        required=True,
        default='scatterplot-sample_swap',
        help='CSV file.'
    )

    options = parser.parse_args()

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    # Read Chr X and Chr Y genes
    X = pd.read_csv(options.X)
    Y = pd.read_csv(options.Y)
    adata.var['X_chr-gene'] = np.in1d(adata.var.index, X['Gene stable ID'])
    adata.var['Y_chr-gene'] = np.in1d(adata.var.index, Y['Gene stable ID'])
    adata.obs['X_chr-sum'] = adata[:, adata.var['X_chr-gene']].X.todense(
    ).sum(axis=1)
    adata.obs['Y_chr-sum'] = adata[:, adata.var['Y_chr-gene']].X.todense(
    ).sum(axis=1)
    df = adata.obs[['sanger_sample_id', 'sex', 'Y_chr-sum', 'X_chr-sum']]
    df = df.groupby(['sanger_sample_id', 'sex']).mean().dropna().reset_index()

    # Save scatterplot with mean expression per sample
    plt = plt9.ggplot(df) + plt9.aes(x='X_chr-sum', y='Y_chr-sum', color='sex')
    plt = plt + plt9.geom_point()
    plt = plt + plt9.ylab("Mean Y chr gene expression ")
    plt = plt + plt9.xlab("Mean X chr gene expression ")
    plt9.ggsave(plt, filename='{}.png'.format(options.o))


if __name__ == '__main__':
    main()
