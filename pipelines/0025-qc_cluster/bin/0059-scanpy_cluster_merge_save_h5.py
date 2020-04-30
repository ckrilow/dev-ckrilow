#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-04-29'
__version__ = '0.0.1'


import argparse
import pandas as pd
import scanpy as sc
import h5py


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object with cluster information. Saves normalized \
            matrix and cluster information as h5 and csv files.
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
        '-od', '--output_dir',
        action='store',
        dest='od',
        default='tenx_from_adata',
        help='Basename of output directory.\
            (default: %(default)s)'
    )

    opt = parser.parse_args()

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=opt.h5)

    df = pd.DataFrame(adata._X)
    df.index = adata.obs.index
    df.columns = adata.var.index

    # Save normalized matrix
    out_f = '{}_X.h5'.format(opt.od)
    hf = h5py.File(out_f, 'w')
    hf.create_dataset('df', data=df, compression='gzip')
    hf.create_dataset('genes', data=df.columns, compression='gzip')
    hf.create_dataset('cells', data=df.index, compression='gzip')
    hf.close()

    # Save metadata with cluster information
    out_f1 = '{}obs.csv'.format(opt.od)
    df1 = pd.DataFrame(adata.obs)
    df1.to_csv(out_f1, sep='\t')


if __name__ == '__main__':
    main()
