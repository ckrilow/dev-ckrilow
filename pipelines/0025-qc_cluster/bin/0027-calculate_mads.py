#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-06-29'
__version__ = '0.0.1'

import argparse
from scipy import stats
import pandas as pd
import plotnine as plt9
import numpy as np
import scanpy as sc


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Calculate one, two or three MADs above median for qc metric \
            specified in key. Saves tsv file with calculated MADs and \
            distribution plot of qc key metric.
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
        '-qck', '--qc_key',
        action='store',
        dest='qc_key',
        default='pct_counts_mito_gene',
        help='QC key for MADs calculation.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output png file. Will have -mads-<qc_key> appended\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    if options.of == '':
        output_file = 'mads-{}'.format(
            options.qc_key
        )
    else:
        output_file = '{}-mads-{}'.format(
            options.of,
            options.qc_key
        )

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    df = adata.obs[[options.qc_key]]
    med = np.median(df[options.qc_key])
    mad = stats.median_abs_deviation(df[options.qc_key])
    mad1 = med+(1*mad)
    mad2 = med+(2*mad)
    mad3 = med+(3*mad)

    dat = {
        'stats': [med, mad1, mad2, mad3],
        'round': [round(med), round(mad1), round(mad2), round(mad3)]}
    df1 = pd.DataFrame(
        dat,
        index=['median', 'median + 1*mad', 'median + 2*mad', 'median + 3*mad'])

    df1.to_csv(
        '{}-mads.tsv'.format(output_file),
        sep='\t',
        index=True
    )

    gplt = plt9.ggplot(df, plt9.aes(x=options.qc_key))
    gplt = gplt + plt9.geom_histogram()
    gplt = gplt + plt9.geom_vline(
                            xintercept=med, linetype="dashed", color="red")
    gplt = gplt + plt9.geom_vline(
                            xintercept=mad1, linetype="dashed", color="red")
    gplt = gplt + plt9.geom_vline(
                            xintercept=mad2, linetype="dashed", color="red")
    gplt = gplt + plt9.geom_vline(
                            xintercept=mad3, linetype="dashed", color="red")

    gplt.save(
        '{}.png'.format(output_file),
        dpi=300,
        width=5,
        height=4
    )


if __name__ == '__main__':
    main()
