#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-05-01'
__version__ = '0.0.1'

import argparse
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
import plotnine as plt9
import csv


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Fits logistic regression to predict labels in adata.obs["cluster"]'
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-tr', '--test_result_tsv',
        action='store',
        dest='test_result',
        required=True,
        help='List of tab-delimited files of logistic regression test_results\
            List should be split by "::" (e.g. file1.tsv.gz::file2.tsv.gz).'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='logistic_regression_resolution_auc',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get the out file base.
    out_file_base = options.of

    # Read in a list of model_report.tsv.gz files.
    # NOTE: cluster id will not necissarily be comparable across reports.
    files = options.model_reports.split('::')

    list_modelreport = []
    for i in range(len(files)):
        _tmp = pd.read_csv(
            files[i],
            sep='\t'
        )
        _tmp['file'] = files[i]
        list_modelreport.append(_tmp)

    # Make one long dataframe.
    df_modelreport = pd.concat(list_modelreport)

    # Save the results
    df_modelreport.to_csv(
        '{}.tsv.gz'.format(out_file_base),
        sep='\t',
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression='gzip'
    )

    # TODO: check for multiple reduced dims or something.

    # Plot the data.
    len_x = len(np.unique(df_modelreport['resolution']))
    len_x2 = len(np.unique(df_modelreport['sparsity']))
    gplt = plt9.ggplot(df_modelreport, plt9.aes(
        fill='sparsity',
        x='resolution',
        y='AUC',
    ))
    gplt = gplt + plt9.theme_bw(base_size=12)
    gplt = gplt + plt9.geom_boxplot(alpha=0.9)
    gplt = gplt + plt9.geom_jitter(
        plt9.aes(color='sparsity'),
        alpha=0.45
    )
    # gplt = gplt + plt9.facet_grid('{} ~ .'.format(label))
    gplt = gplt + plt9.labs(
        x='Resolution',
        y='AUC',
        title=''
    )
    gplt = gplt + plt9.theme(
        # legend_position='none',
        axis_text_x=plt9.element_text(angle=-45, hjust=0)
    )
    if len_x2 != 0 and len_x2 < 9:
        gplt = gplt + plt9.scale_fill_brewer(
            palette='Dark2',
            type='qual'
        )
    gplt.save(
        '{}-resolution_auc.pdf'.format(out_file_base),
        dpi=300,
        width=4*((len_x+len_x2)/4),
        height=5,
        limitsize=False
    )


if __name__ == '__main__':
    main()
