#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import pandas as pd
import csv


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Subsets PCA file for analysis.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-pc', '--tsv_pcs',
        action='store',
        dest='pc',
        default='',
        help='Tab-delimited file of PCs for each cell. First column is\
            cell_barcode. Subsequent columns are PCs.'
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
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <tsv_pcs>-subset)'
    )

    options = parser.parse_args()

    # Fixed settings.
    verbose = True

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-subset'.format(
            os.path.basename(options.pc.rstrip('.tsv.gz'))
        )

    # Load the PCs.
    df_pca = pd.read_csv(options.pc, sep='\t', index_col='cell_barcode')

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

    # Subset down to these PCs.
    df_pca = df_pca.iloc[:, range(0, n_pcs)]

    # Save the clustered data to a data frame.
    df_pca.to_csv(
        '{}.tsv.gz'.format(out_file_base),
        sep='\t',
        index=True,
        quoting=csv.QUOTE_NONNUMERIC,
        index_label='cell_barcode',
        na_rep='',
        compression='gzip'
    )


if __name__ == '__main__':
    main()
