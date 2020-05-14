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


def _make_plots(
    df_plt,
    out_file_base,
    facet_grid=''
):
    # Plot AUC.
    len_x = len(np.unique(df_plt['resolution']))
    if 'sparsity_l1' in df_plt.columns:
        df_plt['Sparsity'] = df_plt['sparsity_l1']
        len_x2 = len(np.unique(df_plt['Sparsity']))
    else:
        len_x2 = 0
    if len_x2 > 1:
        gplt = plt9.ggplot(df_plt, plt9.aes(
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
        gplt = plt9.ggplot(df_plt, plt9.aes(
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
    if facet_grid != '':
        gplt = gplt + plt9.facet_grid('{} ~ .'.format(facet_grid))
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
        gplt = plt9.ggplot(df_plt, plt9.aes(
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
        gplt = plt9.ggplot(df_plt, plt9.aes(
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
    if facet_grid != '':
        gplt = gplt + plt9.facet_grid('{} ~ .'.format(facet_grid))
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

    # get all pairwise combinations of resolutions and neighbors
    resolutions = df_modelreport['resolution'].unique()
    n_neighbors = df_modelreport['neighbors__n_neighbors'].unique()
    iters = [(x, y) for x in resolutions for y in n_neighbors]

    # Check the data.
    # Also add flag if mulitple neighbors have been used for clustering
    not_unique_columns = [
        'cell_label', 'precision', 'recall', 'f1-score',
        'support', 'AUC', 'n_cells_full_dataset', 'n_cells_training_dataset',
        'is_cluster'
    ]
    unique_columns = [
        i for i in df_modelreport.columns if i not in not_unique_columns
    ]
    for i__res, i__n_neighbors in iters:
        # Get the subset of data
        df_tmp = df_modelreport[df_modelreport['resolution'] == i__res]
        df_tmp = df_tmp[df_tmp['neighbors__n_neighbors'] == i__n_neighbors]
        # Check to make sure we have the right data ... for every non-class
        # column, there should be only one value. Otherwise, we are not
        # plotting a unique run.
        for index, value in df_tmp[unique_columns].nunique().items():
            if 'cell_label' not in index and value > 1:
                print(index, value)
                raise Exception('ERROR: multiple resolution files combined.')

    # Make sure we just get cluster variables as there is other info included
    # in model reports.
    df_modelreport = df_modelreport[df_modelreport['is_cluster']]

    # If we have multiple neighbor values, make a facet plot.
    if len(n_neighbors) > 1:
        _make_plots(
            df_modelreport,
            out_file_base,
            facet_grid='neighbors__n_neighbors'
        )
    # Also make seperate plots for each n_neighbors value.
    for i__n_neighbors in df_modelreport['neighbors__n_neighbors'].unique():
        df_tmp = df_modelreport[
            df_modelreport['neighbors__n_neighbors'] == i__n_neighbors
        ]
        _make_plots(
            df_tmp,
            '{}-n_neighbors={}'.format(out_file_base, i__n_neighbors),
            facet_grid=''
        )


if __name__ == '__main__':
    main()
