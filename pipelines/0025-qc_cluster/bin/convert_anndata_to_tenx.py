#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
from distutils.version import LooseVersion
import scipy.io as sio
import gzip
import pandas as pd
import scanpy as sc


def anndata_to_tenx(
    adata,
    raw=False,
    out_file='',
    out_dir='tenx_from_adata'
):
    """Write 10x like data from anndata.

    Parameters
    ----------
    adata : pandas.DataFrame
        Description of parameter `adata`.
    raw : boolean
        Description of parameter `raw`.
    output_file : string
        Description of parameter `output_file`.
    out_dir : string
        Description of parameter `out_dir`.

    Returns
    -------
    execution_code : int
    """
    # Make the output directory if it does not exist.
    if out_dir == '':
        out_dir = os.getcwd()
    else:
        os.makedirs(out_dir)

    # If raw, overwrite adata to be the raw dataframe.
    if raw:
        adata = adata.raw.to_adata()

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # Save the barcodes.
    pd.DataFrame(adata.obs.index).to_csv(
        os.path.join(
            out_dir,
            '{}-barcodes.tsv.gz'.format(out_file)
        ),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )

    # Save the features.
    for i in ['gene_symbols', 'feature_types']:
        if i not in adata.var.columns:
            raise Exception('Could not find expected columns in adata.var.')
    df_features = adata.var.loc[:, ['gene_symbols', 'feature_types']]
    df_features.to_csv(
        os.path.join(
            out_dir,
            '{}-features.tsv.gz'.format(out_file)
        ),
        sep='\t',
        compression=compression_opts,
        index=True,
        header=False
    )

    # Save the metadata.
    adata.obs.to_csv(
        os.path.join(
            out_dir,
            '{}-metadata.tsv.gz'.format(out_file)
        ),
        sep='\t',
        compression=compression_opts,
        index=True,
        header=True,
        index_label='cell_barcode'
    )

    # Save all matricies - both the core matrix (X) as well as any layers.
    matrices = set(['X']).union(set(adata.layers.keys()))
    for mtx_i in matrices:
        if mtx_i == 'X':
            out_mtx = adata.X.transpose()
        else:
            out_mtx = adata.layers[mtx_i].transpose()
        out = os.path.join(
            out_dir,
            '{}-matrix-{}.mtx.gz'.format(out_file, mtx_i)
        )
        with gzip.open(out, 'wb', compresslevel=9) as f:
            sio.mmwrite(
                f,
                out_mtx,
                comment='metadata_json: {"format_version": 2}',
                field='real'  # can be integer if counts otherwise real TODO
            )

    return 0


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read an H5 AnnData object and write matrix files similar to 10x
            output. The results can be read into R via Seurat::Read10X and
            Seurat::CreateSeuratObject. This script requires pandas >1.0.1
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
        '-r', '--raw',
        action='store_true',
        dest='r',
        default=False,
        help='Convert data in the raw anndata slot. Note that data in the\
            raw slot may have different dimensions than data in the main\
            slot. Therefore, ensure that the output from this script run\
            with the --raw flag does not go to the same directory as it\
            called without the --raw flag since the barcodes and features\
            files may not properly match in the output.\
            (default: %(default)s)'
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
    adata = sc.read_h5ad(filename=options.h5)

    # Run the conversion function.
    _ = anndata_to_tenx(
        adata,
        raw=options.r,
        out_file=options.of,
        out_dir=options.od
    )


if __name__ == '__main__':
    main()
