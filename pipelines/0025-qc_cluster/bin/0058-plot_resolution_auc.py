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
            Plots AUC from model_report.tsv.gz
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-mr', '--model_reports',
        action='store',
        dest='model_reports',
        required=True,
        help='List of tab-delimited files of model_reports.\
            List should be split by "::" (e.g. file1.tsv.gz::file2.tsv.gz).'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='resolution_auc',
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
        '{}-merged_model_report.tsv.gz'.format(out_file_base),
        sep='\t',
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression='gzip'
    )

    if 'resolution' not in df_modelreport.columns:
        if 'cluster__resolution' in df_modelreport.columns:
            df_modelreport['resolution'] = df_modelreport[
                'cluster__resolution'
            ]
        # else:
        #     df_modelreport['resolution'] = df_modelreport['file']

    # Make sure resolution is set to a categorical variable.
    df_modelreport['resolution'] = df_modelreport['resolution'].astype(
        'category'
    )

    # Make sure we just get cluster variables as there is other info included
    # in model reports.
    df_modelreport = df_modelreport[df_modelreport['is_cluster']]

    # Plot AUC.
    len_x = len(np.unique(df_modelreport['resolution']))
    if 'sparsity_l1' in df_modelreport.columns:
        df_modelreport['Sparsity'] = df_modelreport['sparsity_l1']
        len_x2 = len(np.unique(df_modelreport['Sparsity']))
    else:
        len_x2 = 0
    if len_x2 > 1:
        gplt = plt9.ggplot(df_modelreport, plt9.aes(
            fill='Sparsity',
            x='resolution',
            y='AUC',
        ))
        gplt = gplt + plt9.geom_boxplot(
            alpha=0.8,
            outlier_alpha=0
        )
        gplt = gplt + plt9.geom_jitter(
            plt9.aes(color='Sparsity'),
            alpha=0.25,
            width=0.2
        )
    else:
        gplt = plt9.ggplot(df_modelreport, plt9.aes(
            x='resolution',
            y='AUC',
        ))
        gplt = gplt + plt9.geom_boxplot(
            alpha=0.8,
            outlier_alpha=0
        )
        gplt = gplt + plt9.geom_jitter(
            alpha=0.25,
            width=0.2
        )
    gplt = gplt + plt9.theme_bw(base_size=12)
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

    # Plot F1
    if len_x2 > 1:
        gplt = plt9.ggplot(df_modelreport, plt9.aes(
            fill='Sparsity',
            x='resolution',
            y='f1-score',
        ))
        gplt = gplt + plt9.geom_boxplot(
            alpha=0.8,
            outlier_alpha=0
        )
        gplt = gplt + plt9.geom_jitter(
            plt9.aes(color='Sparsity'),
            alpha=0.25,
            width=0.2
        )
    else:
        gplt = plt9.ggplot(df_modelreport, plt9.aes(
            x='resolution',
            y='f1-score',
        ))
        gplt = gplt + plt9.geom_boxplot(
            alpha=0.8,
            outlier_alpha=0
        )
        gplt = gplt + plt9.geom_jitter(
            alpha=0.25,
            width=0.2
        )
    gplt = gplt + plt9.theme_bw(base_size=12)
    # gplt = gplt + plt9.facet_grid('{} ~ .'.format(label))
    gplt = gplt + plt9.labs(
        x='Resolution',
        y='F1 score',
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
        '{}-resolution_f1score.pdf'.format(out_file_base),
        dpi=300,
        width=4*((len_x+len_x2)/4),
        height=5,
        limitsize=False
    )


if __name__ == '__main__':
    main()
