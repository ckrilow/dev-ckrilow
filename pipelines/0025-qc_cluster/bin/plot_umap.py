#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# Silence NumbaPerformanceWarning in umap. See below:
# https://github.com/lmcinnes/umap/issues/252
import warnings
from numba.errors import NumbaPerformanceWarning
warnings.filterwarnings('ignore', category=NumbaPerformanceWarning)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and PCs file. Generates UMAP.
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
        '-cq', '--colors_quantitative',
        action='store',
        dest='cq',
        default='',
        help='Comma seperated list of quantitative variable names for colors.\
            (default: "")'
    )

    parser.add_argument(
        '-cc', '--colors_categorical',
        action='store',
        dest='cc',
        default='',
        help='Comma seperated list of categorical variable names for colors.\
            (default: "")'
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
        '-nn', '--n_neighbors',
        action='store',
        dest='n_neighbors',
        default=15,
        type=int,
        help='Number of neighbors for sc.pp.neighbors call\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-umd', '--umap_min_dist',
        action='store',
        dest='umap_min_dist',
        default=0.5,
        type=float,
        help='The effective minimum distance between embedded points. Smaller\
            values will result in a more clustered/clumped embedding where\
            nearby points on the manifold are drawn closer together, while\
            larger values will result on a more even dispersal of points.\
            The value should be set relative to the spread value, which\
            determines the scale at which embedded points will be spread out.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-us', '--umap_spread',
        action='store',
        dest='umap_spread',
        default=1.0,
        type=float,
        help='The effective scale of embedded points. In combination with\
            min_dist this determines how clustered/clumped the embedded\
            points are.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-uip', '--umap_init_pos',
        action='store',
        dest='umap_init_pos',
        default='X_pca',
        help='How to initialize the low dimensional embedding.\
            Valid options: any key for adata.obsm,\
            ’paga’: positions from paga(),\
            ’spectral’: use a spectral embedding of the graph,\
            ’random’: assign initial embedding positions at random.\
            (default: X_pca, the slot where tsv_pcs is stored if provided)'
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
            (default: <h5_anndata>-<tsv_pcs>-umap)'
    )

    options = parser.parse_args()

    # Fixed settings.
    verbose = True

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    sc.set_figure_params(dpi_save=300)

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
    # df_pca = pd.read_csv(
    #     'adata-pcs-harmony.tsv.gz',
    #     sep='\t',
    #     index_col='cell_barcode'
    # )

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
    # Subset number of PCs to be exactly nPCs - here we assume PCs are ordered.
    print('Subetting PCs - we assume they are ordered by column index')
    df_pca = df_pca.iloc[:, range(0, n_pcs)]
    print('PC columns:\t{}'.format(np.array_str(df_pca.columns)))

    # Add the reduced dimensions to the AnnData object.
    adata.obsm['X_pca'] = df_pca.loc[adata.obs.index, :].values.copy()

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-{}-umap'.format(
            os.path.basename(options.h5.rstrip('.h5ad')),
            os.path.basename(options.pc.rstrip('.tsv.gz'))
        )
    # Append the parameters to the output file.
    out_file_base = '{}.number_pcs={}'.format(
        out_file_base,
        n_pcs
    )
    out_file_base = '{}.n_neighbors={}'.format(
        out_file_base,
        options.n_neighbors
    )
    out_file_base = '{}.umap_min_dist={}'.format(
        out_file_base,
        str(options.umap_min_dist).replace('.', 'pt')
    )
    out_file_base = '{}.umap_spread={}'.format(
        out_file_base,
        str(options.umap_spread).replace('.', 'pt')
    )
    out_file_base = '{}.umap_init_pos={}'.format(
        out_file_base,
        options.umap_init_pos
    )

    # Parse the color variables.
    colors_quantitative = []
    if options.cq != '':
        colors_quantitative = options.cq.split(',')

    colors_categorical = []
    if options.cc != '':
        colors_categorical = options.cc.split(',')

    if len(colors_quantitative) == 0 and len(colors_categorical) == 0:
        raise Exception('Specify a color value.')

    # Nice large palette.
    colors_large_palette = [
        '#0F4A9C', '#3F84AA', '#C9EBFB', '#8DB5CE', '#C594BF', '#DFCDE4',
        '#B51D8D', '#6f347a', '#683612', '#B3793B', '#357A6F', '#989898',
        '#CE778D', '#7F6874', '#E09D37', '#FACB12', '#2B6823', '#A0CC47',
        '#77783C', '#EF4E22', '#AF1F26'
    ]
    # Add colors_large_palette to adata.uns.
    # adata.uns["annotation_colors"] = colors_large_palette

    # Calculate neighbors for on the specified PCs.
    sc.pp.neighbors(
        adata,
        use_rep='X_pca',
        n_pcs=n_pcs,
        n_neighbors=options.n_neighbors,  # Scanpy default = 15
        copy=False
    )

    # If init with paga, plot paga first - NOTE we can only do this if
    if options.umap_init_pos == 'paga' and 'paga' not in adata.uns:
        print(
            'Trying to call sc.tl.paga.',
            'NOTE: requires one to have clustered the data.'
        )
        sc.tl.paga(
            adata,
            use_rna_velocity=False,
            copy=False
        )

    # UMAP
    sc.tl.umap(
        adata,
        min_dist=options.umap_min_dist,  # Scanpy default = 0.05
        spread=options.umap_spread,  # Scanpy default = 1.0
        init_pos=options.umap_init_pos,  # Scanpy default = spectral
        copy=False
    )

    # NOTE: If the color var is a gene, you should color by ln(CPM+1).
    #       By default these sc.pl.umap uses the .raw attribute of AnnData
    #       if present which is assumed to be ln(CPM+1).

    # For each variable, loop over and set color accordingly. Save
    # the results.
    for var in colors_quantitative:
        fig = sc.pl.umap(
            adata,
            color=var,
            alpha=0.25,
            return_fig=True
        )
        fig.savefig(
            '{}-{}.png'.format(out_file_base, var),
            dpi=300,
            bbox_inches='tight'
        )
    for var in colors_categorical:
        print(var)
        # Cast to category - required for booleans.
        adata.obs[var] = adata.obs[var].astype('category')
        n_categories = len(adata.obs[var].cat.categories)
        color_palette = None
        if n_categories <= len(plt.get_cmap('Dark2').colors):
            color_palette = 'Dark2'
        elif n_categories <= len(colors_large_palette):
            color_palette = colors_large_palette
        fig = sc.pl.umap(
            adata,
            color=var,
            palette=color_palette,
            alpha=0.25,
            return_fig=True
        )
        fig.savefig(
            '{}-{}.png'.format(out_file_base, var),
            dpi=300,
            bbox_inches='tight'
        )

    # In some ocassions, you might still observe disconnected clusters and
    # similar connectivity violations. They can usually be remedied by running:
    # sc.tl.paga(adata)
    # From below, remove `plot=False` if you want to see the coarse-grained
    # graph
    # sc.pl.paga(adata, plot=False)
    # sc.tl.umap(adata, init_pos='paga')


if __name__ == '__main__':
    main()
