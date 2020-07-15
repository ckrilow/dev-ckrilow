#!/usr/bin/env python


__author__ = 'Henry Taylor'
__date__ = '2020-07-09'
__version__ = '0.0.1'

import argparse
import os
from distutils.version import LooseVersion
import scipy
import scipy.io
import gzip
import pandas as pd
import scanpy as sc


def cellbender_to_tenxmatrix(
    adata,
    out_file='',
    out_dir='tenx_from_adata',
    verbose=True
):
    """Write 10x like data from 10x H5.

    Parameters
    ----------
    adata : pandas.DataFrame
        Description of parameter `adata`.
    output_file : string
        Description of parameter `output_file`.
    out_dir : string
        Description of parameter `out_dir`.
    verbose : boolean
        Description of parameter `verbose`.

    Returns
    -------
    execution_code : int
    """
    # Make the output directory if it does not exist.
    if out_dir == '':
        out_dir = os.getcwd()
    else:
        os.makedirs(out_dir, exist_ok=True)

    # Set up out_file
    if out_file != '':
        out_file = '{}-'.format(out_file)

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Save the barcodes.
    out_f = os.path.join(
        out_dir,
        '{}barcodes.tsv.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))
    pd.DataFrame(adata.obs.index).to_csv(
        out_f,
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )

    # Save the features.
    # CellBender output: gene_symbols as index, enseble as `gene_ids`
    # Type is stored in adata.var.type
    out_f = os.path.join(
        out_dir,
        '{}features.tsv.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))

    # Check if read in by gene_symbol or gene_id
    if 'gene_symbol' in adata.var.columns:
        df_features = pd.DataFrame(data=[ adata.var.index.values, adata.var.loc[:,'gene_symbol'].values, ["Gene Expression"] * adata.n_vars ]).T
    elif 'gene_ids' in adata.var.columns:
        df_features = pd.DataFrame(data=[ adata.var.loc[:,'gene_ids'].values, adata.var.index.values, ["Gene Expression"] * adata.n_vars ]).T
    else:
        raise Exception('Could not find "gene_symbols" or "gene_ids" in adata.var')

    df_features.to_csv(
        out_f,
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )

    # Save the count-adjusted matrix
    out_mtx = adata.X.transpose()
    if not isinstance(out_mtx, scipy.sparse.csr.csr_matrix):
        out_mtx = scipy.sparse.csr_matrix(out_mtx)
    out_f = os.path.join(
        out_dir,
        '{}matrix.mtx.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))
        # print(out_mtx.dtype)  # it looks like even counts stored as float
    with gzip.open(out_f, 'wb', compresslevel=9) as f:
        scipy.io.mmwrite(
            f,
            out_mtx,
            comment='metadata_json: {"format_version": 2}'
            # field='real'  # can be integer if counts otherwise real
        )

    if verbose:
        print(
            'Done.'
        )

    return 0


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read an H5 10x-formatted object and write matrix files similar to 10x
            output. This script requires pandas >1.0.1
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5', '--h5_10x',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-g', '--genome',
        action='store',
        dest='g',
        default="background_removed",
        help='Filter expression to genes within this genome. For CellBender,\
        this will be "background_removed"'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files.\
            (default: no tag basename)'
    )

    parser.add_argument(
        '-od', '--output_dir',
        action='store',
        dest='od',
        default='tenx_from_adata',
        help='Basename of output directory.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Load the AnnData file
    adata = sc.read_10x_h5(filename=options.h5, genome=options.g)

    # Run the conversion function.
    _ = cellbender_to_tenxmatrix(
        adata,
        out_file=options.of,
        out_dir=options.od
    )


if __name__ == '__main__':
    main()
