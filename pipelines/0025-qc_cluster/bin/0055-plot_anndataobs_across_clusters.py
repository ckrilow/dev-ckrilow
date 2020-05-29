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
        default='',
        help='Pheno column to be boxplotted by cluster.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='plot_boxplot_cluster',
        help='Basename of output png file. Will have .png appended.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    adata = sc.read_h5ad(filename=options.h5)

    pheno_to_plot = options.pheno_columns.split(',')

    # Plot the data.
    for pheno in pheno_to_plot:
        df = adata.obs[[pheno]].copy()
        df['cluster'] = adata.obs['cluster'].copy()
        gplt = gplt9.ggplot(df)
        gplt = gplt + gplt9.geom_boxplot(gplt9.aes(x='cluster', y=pheno))
        gplt = gplt + gplt9.theme(axis_text_x=gplt9.element_text(angle=90))
        gplt.save(
            'boxplot-{}.png'.format(pheno),
            dpi=300
        )
        # Add log10 transformation
        gplt = gplt + gplt9.scale_y_continuous(
            # trans='log10',
            # labels=comma_labels,
            minor_breaks=0
        )
        gplt.save(
            'boxplot_log10-{}.png'.format(pheno),
            dpi=300
        )


if __name__ == '__main__':
    main()
