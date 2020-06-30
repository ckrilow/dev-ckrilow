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

    output_file = options.of
    if output_file == '':
        output_file = 'mads-{}'.format(
            options.qc_key
        )

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    qc_keys = options.qc_key.split(',')
    df_dict = {
        'variable': [],
        'experiment_id': [],
        'cutoff_type': [],
        'cutoff': [],
        'cutoff_round': []
    }
    for qc_key in qc_keys:
        med = np.median(adata.obs[qc_key])
        mad = stats.median_abs_deviation(adata.obs[qc_key])
        for i in [1, 2, 3]:
            df_dict['variable'].append(qc_key)
            df_dict['experiment_id'].append('all_experiment_ids')
            df_dict['cutoff_type'].append('median+({}*mad)'.format(i))
            df_dict['cutoff'].append(med+(i*mad))
            df_dict['cutoff_round'].append(round(med+(i*mad)))
        for i in [1, 2, 3]:
            df_dict['variable'].append(qc_key)
            df_dict['experiment_id'].append('all_experiment_ids')
            df_dict['cutoff_type'].append('median-({}*mad)'.format(i))
            df_dict['cutoff'].append(med-(i*mad))
            df_dict['cutoff_round'].append(round(med-(i*mad)))

    df = pd.DataFrame(df_dict)
    df.to_csv(
        '{}.tsv'.format(output_file),
        sep='\t',
        index=False
    )

    for qc_key in qc_keys:
        gplt = plt9.ggplot(adata.obs, plt9.aes(x=qc_key))
        gplt = gplt + plt9.geom_histogram()
        filt = (df['variable'] == qc_key) & np.in1d(df['cutoff_type'].values, [
            'median+({}*mad)'.format(i) for i in [1, 2, 3]
        ])
        for cut in df.loc[filt, 'cutoff']:
            gplt = gplt + plt9.geom_vline(
                xintercept=cut,
                linetype="dashed",
                color="red"
            )
        gplt.save(
            '{}-{}.png'.format(output_file, qc_key),
            dpi=300,
            width=5,
            height=4
        )


if __name__ == '__main__':
    main()
