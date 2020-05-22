#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import warnings

# Silence NumbaPerformanceWarning in umap. See below:
# https://github.com/lmcinnes/umap/issues/252
from numba.errors import NumbaPerformanceWarning
warnings.filterwarnings('ignore', category=NumbaPerformanceWarning)

# To error out when 'FloatingPointError: divide by zero encountered in power'
# from UMAP. This is usually caused by a poor choice of min_dist and spread
# parameters.
np.seterr(all='raise')

# NOTE: There are two ways we could write this script.
# Version 1:
# 1. Save pickle of connectivities dict.
# 2. Save umap params tsv files.
# 3. Save umap matrix to tsv.
#
# Version 2:
# 1. Save adata
#
# We are doing #2 because the merge of many UMAPS saved in adata will be
# easier.


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and PCs file. Calculates UMAP given parameters.
            Saves UMAP object.
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
        help='Number of neighbors for sc.pp.neighbors call.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-uinit', '--umap_init',
        action='store',
        dest='umap_init',
        default='X_pca',
        help='How to initialize the low dimensional embedding.\
            Valid options: any key for adata.obsm,\
            ’paga’: positions from paga(),\
            ’spectral’: use a spectral embedding of the graph,\
            ’random’: assign initial embedding positions at random.\
            (default: X_pca, the slot where tsv_pcs is stored if provided)'
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
        help='The minimum distance apart that points are allowed to be in the\
            low dimensional representation (effective scale of embedded points\
            ). In combination with min_dist this determines how\
            clustered/clumped the embedded points are.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-cd', '--calculate_densities',
        action='store_true',
        dest='calculate_densities',
        default=False,
        help='Calculate point density.\
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
        '--force_recalculate_neighbors',
        action='store_true',
        dest='calculate_neighbors',
        default=False,
        help='Calculate neighbor graph even if it already exists in the\
            AnnData (which it my do so if you already ran BBKNN).\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <h5_anndata>-<tsv_pcs>.h5ad)'
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
    # NOTE: In the BKNN instnace, we could totally ignore this, since the
    #       neighbors graph is already calculated.
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

    # Add the reduced dimensions to the AnnData object.
    # NOTE: We do this later... only if we need to.
    # adata.obsm['X_pca'] = df_pca.loc[adata.obs.index, :].values.copy()

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
    print('Subetting PCs - we assume they are ordered by column index.')
    df_pca = df_pca.iloc[:, range(0, n_pcs)]
    print('PC columns:\t{}'.format(np.array_str(df_pca.columns)))

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-{}-umap'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.')),
            os.path.basename(options.pc.rstrip('tsv.gz').rstript('.'))
        )

    # Parse the neighbors iterations.
    i__n_neighbors = options.n_neighbors
    # Parse the min_dist iterations.
    i__min_dist = options.umap_min_dist
    # Parse the neighbors iterations.
    i__spread = options.umap_spread

    # Check input parameters
    if not (2 <= i__n_neighbors <= 100):
        # Recommended in parameter documentation:
        # https://umap-learn.readthedocs.io/en/latest/api.html
        warnings.warn(
            'WARNING: it is suggested to set n_neighbors to a {}'.format(
                'value between 2-100.'
            )
        )
    if not (0.0 <= i__min_dist <= 1.0):
        # Recommended here: https://github.com/lmcinnes/umap/issues/249
        warnings.warn(
            'WARNING: it is suggested to set umap_min_dist to a {}'.format(
                'value between 0-1.'
            )
        )
    if not (0.0 <= i__spread <= 3.0):
        # Recommendation based on single cell experience.
        warnings.warn(
            'WARNING: it is suggested to set umap_spread to a {}'.format(
                'value between 0-3.'
            )
        )

    # Calculate neighbors for on the specified PCs.
    # By default saved to adata.uns['neighbors']
    #
    # First, however, check to see if adata.uns['neighbors'] already exists
    # ...and unless the user tells us not to, use that slot, not calculating
    # neighbors. This default behaviour is to accommodate the instance when
    # bbknn has been run on the data.
    if 'neighbors' not in adata.uns or options.calculate_neighbors:
        # Add the reduced dimensions to the AnnData object.
        adata.obsm['X_pca'] = df_pca.loc[adata.obs.index, :].values.copy()

        sc.pp.neighbors(
            adata,
            use_rep='X_pca',
            n_pcs=n_pcs,
            n_neighbors=i__n_neighbors,  # Scanpy default = 15
            copy=False
        )
    else:
        warnings.warn(
            'WARNING: found neighbors slot in adata.uns. {}'.format(
                'Not calculating neighbors (ignoring n_neighbors parameter).'
            )
        )
        # If we are using the pre-calculated neighbors drop npcs note.
        # if 'n_pcs' in adata.uns['neighbors']['params']:
        n_pcs = adata.uns['neighbors']['params']['n_pcs']
        i__n_neighbors = adata.uns['neighbors']['params']['n_neighbors']
    # adata.uns['neighbors__{}'.format(plt__label)] = adata.uns['neighbors']

    # TODO: add paga
    # # If init with paga, plot paga first - NOTE we can only do this if
    # if options.umap_init == 'paga' and 'paga' not in adata.uns:
    #     print(
    #         'Trying to call sc.tl.paga.',
    #         'NOTE: requires one to have clustered the data.'
    #     )
    #     sc.tl.paga(
    #         adata,
    #         use_rna_velocity=False,
    #         copy=False
    #     )

    # UMAP
    # Saved to adata.uns['umap'] and adata.obsm['X_umap']
    sc.tl.umap(
        adata,
        min_dist=i__min_dist,  # Scanpy default = 0.05
        spread=i__spread,  # Scanpy default = 1.0
        init_pos=options.umap_init,  # Scanpy default = spectral
        # For some reason cannot access neighbors key slot, thus we
        # must keep uns['neighbors'] until we have run this.
        # neighbors_key='neighbors__{}'.format(plt__label),
        copy=False
    )
    # Add umap info to params stash
    # adata.uns['umap']['params']['tsv_reduced_dims'] = options.pc
    adata.uns['umap']['params']['n_reduced_dims_input'] = n_pcs
    adata.uns['umap']['params']['n_neighbors'] = i__n_neighbors
    adata.uns['umap']['params']['umap_init'] = options.umap_init
    adata.uns['umap']['params']['umap_min_dist'] = i__min_dist
    adata.uns['umap']['params']['umap_spread'] = i__spread
    # adata.uns['umap']['params']['density'] = ''

    if options.calculate_densities:
        sc.tl.embedding_density(
            adata,
            basis='umap'
        )
        # Density is the root string of the density parameters.
        # If no density calculations, then this will be ''
        adata.uns['umap']['params']['density'] = 'umap_density'
        # Rename density estimates
        # adata.obs[
        #     'umap__{}__density'.format(plt__label)
        # ] = adata.obs.pop('umap_density')
        # adata.uns[
        #     'umap__{}__density_params'.format(plt__label)
        # ] = adata.uns.pop('umap_density_params')

    # Rename UMAP params
    # NOTE: we don't do this because scanpy functions by default look for
    # umap params in the adata.uns['umap_params'] slot.
    # adata.uns['umap_params'] = adata.uns.pop('umap')

    adata.write(
        '{}.h5ad'.format(out_file_base),
        compression='gzip'
    )


if __name__ == '__main__':
    main()
