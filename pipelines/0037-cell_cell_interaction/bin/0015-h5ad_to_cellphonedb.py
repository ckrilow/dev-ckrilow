#!/usr/bin/env python


__author__ = 'Julie Matte'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import scanpy as sc
import pandas as pd
import csv


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Makes callphonedb input from h5ad file.
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
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <adata>-cellphonedb)'
    )

    options = parser.parse_args()

    # Fixed settings.
    verbose = True

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-cellphonedb'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.'))
        )

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)


    # If the out file is adata one can save like this.
    # df_pca.to_csv(
    #     '{}.tsv.gz'.format(out_file_base),
    #     sep='\t',
    #     index=True,
    #     quoting=csv.QUOTE_NONNUMERIC,
    #     index_label='cell_barcode',
    #     na_rep='',
    #     compression='gzip'
    # )


if __name__ == '__main__':
    main()
