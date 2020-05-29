#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-05-29'
__version__ = '0.0.1'

import argparse
import os
import numpy as np
import scipy as sp
import scanpy as sc


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and preps for cellxgene pipeline.
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

    # parser.add_argument(
    #     '--force_recalculate_neighbors',
    #     action='store_true',
    #     dest='calculate_neighbors',
    #     default=False,
    #     help='Calculate neighbor graph even if it already exists in the\
    #         AnnData (which it my do so if you already ran BBKNN).\
    #         (default: %(default)s)'
    # )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <h5_anndata>-<tsv_pcs>-cellxgene)'
    )

    options = parser.parse_args()

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-cellxgene'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.'))
        )

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)

    # Format file as described here:
    # https://cellgeni.readthedocs.io/en/latest/visualisations.html

    # Set X to cp10k
    adata.X = np.expm1(adata.layers['log1p_cp10k'])

    # cellxgene works faster when the expression matrix is stored in CSC
    # (compressed sparse column) format instead of CSR
    # (compressed sparse row) format
    adata.X = sp.sparse.csc_matrix(adata.X)

    # Try to set default UMAP using umap_min_dist=1pt0 and umap_spread=1pt0
    umap_default = [i for i in adata.obsm.keys() if 'umap_min_dist=1pt0' in i]
    if len(umap_default) == 0:
        umap_default = [
            i for i in adata.obsm.keys() if 'umap_spread=1pt0' in i
        ]
        if len(umap_default) > 0:
            adata.obsm['X_umap'] = adata.obsm[umap_default.pop()]
    elif len(umap_default) >= 0:
        umap_default2 = [
            i for i in umap_default if 'umap_spread=1pt0' in i
        ]
        if len(umap_default2) > 0:
            adata.obsm['X_umap'] = adata.obsm[umap_default2.pop()]
        else:
            adata.obsm['X_umap'] = adata.obsm[umap_default.pop()]
    # NOTE: an alternative would be to run UMAPs here.

    # For coloring, make sure categorical variables are properly set to
    # Categorical
    # adata.obs['metadata_name'] = pd.Categorical(adata.obs['metadata_name'])
    # adata.obs['metadata_name'] = np.float32(adata.obs['metadata_name'])

    # Save the resulting anndata
    adata.write('{}.h5ad'.format(out_file_base), compression='gzip')


if __name__ == '__main__':
    main()
