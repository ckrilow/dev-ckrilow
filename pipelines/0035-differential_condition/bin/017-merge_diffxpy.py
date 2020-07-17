#!/usr/bin/env python


__author__ = 'Henry Taylor'
__date__ = '2020-07-17'
__version__ = '0.0.1'

import argparse
from distutils.version import LooseVersion
import pandas as pd


def merge_diffxpy_dataframes(
    dataframe_keys,
    dataframe_paths,
    outfile,
    verbose=True
):
    ## Start with first dataframe and get dimensions from df.
    df_merged = pd.read_csv(dataframe_paths[0], sep="\t", header=0) ## Expect a header

    if verbose:
        print("Starting to merge dataframes. Based on the first dataframe "
            "expecting dataframes with {} columns and {} rows.".format(
            len(df_merged.columns),
            len(df_merged)
        ))

    for i in range(1, len(dataframe_keys)):
        key = dataframe_keys[i]
        path = dataframe_paths[i]

        if verbose:
            print("Merging dataframe with the key: '{}'".format(key))

        df = pd.read_csv(path, sep="\t", header=0) ## Expect a header

        len_before_merge = len(df_merged)
        df_merged = df_merged.merge(
            df,
            how="outer",
            left_index=False,
            right_index=False
        )

        if len_before_merge == len(df_merged) and verbose:
            print("Dataframe with the key '{}' is not unique. The length "
                "before merge is the same as the length after merge.".format(key))

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    df_merged.to_csv(
        '{}-de_results.tsv.gz'.format(outfile),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )



def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Merge results of the differential expression tests.
            """
    )

    parser.add_argument(
        '-keys', '--dataframe_keys',
        action='store',
        dest='keys',
        required=True,
        help='Comma separated list of the keys for each dataframe.'
    )

    parser.add_argument(
        '-paths', '--dataframe_paths',
        action='store',
        dest='paths',
        required=True,
        help='Comma separated list of the paths to dataframes. Indices \
            should correspond with the keys.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='diffxpy',
        help='Basename of output file.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get the parameters
    keys = options.keys.split(",")
    paths = options.paths.split(",")

    _ = merge_diffxpy_dataframes(
        keys,
        paths,
        options.of
    )


if __name__ == '__main__':
    main()
