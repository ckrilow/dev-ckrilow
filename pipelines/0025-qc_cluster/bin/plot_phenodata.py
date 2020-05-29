#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-05-28'
__version__ = '0.0.1'

import argparse
import scanpy as sc
import plotnine as plt9


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and list of phenotypes. Plot boxplots of \
            phenotypes across clusters.
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
        '--pheno_columns',
        action='store',
        dest='pheno_columns',
        default='none',
        help='Pheno column to be boxplotted by cluster.'
    )

    options = parser.parse_args()

    adata = sc.read_h5ad(filename=options.h5)

    pheno_to_plot = options.pheno_columns.split(',')
    if len(pheno_to_plot) == 0:
        pheno_to_plot = ['none']
    
    # Plot the data.
    for pheno in pheno_to_plot:
        df = adata.obs[[pheno]].copy()
        df['cluster'] = adata.obs['cluster'].copy()        
        plt = plt9.ggplot(df)
        plt = plt + plt9.geom_boxplot(plt9.aes(x='cluster', y=pheno))
        plt = plt + plt9.theme(axis_text_x=plt9.element_text(angle=90))
        plt.save(
         'boxplot{}.png'.format(pheno),
         dpi=300
         )

if __name__ == '__main__':
    main()
