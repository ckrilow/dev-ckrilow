#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-04-29'
__version__ = '0.0.1'


import argparse
import pandas as pd
import scanpy as sc
import os


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Reads Anndata and clustering progress csv file. Adds clustering \
            progress information to Anndata and saves as h5ad file.
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
        '-m', '--merging_file',
        action='store',
        dest='od',
        help='Merging progress file.'
    )

    opt = parser.parse_args()

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=opt.h5)

    # Bind cluster merging information to Anndata file
    merging_progress = pd.read_table(opt.m)
    merging_progress = merging_progress.apply(lambda x: x.astype('category'))
    adata.obs = pd.concat([adata.obs, merging_progress], axis=1)

    # Save Anndata to .h5ad
    out_file_base = os.path.basename('{}'.format(os.path.splitext(opt.h5)[0]))
    adata.write('{}_merged.h5ad'.format(out_file_base), compression='gzip')

if __name__ == '__main__':
    main()
