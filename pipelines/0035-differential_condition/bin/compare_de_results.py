#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-07-27'
__version__ = '0.0.1'

import argparse
import numpy as np
import pandas as pd
import plotnine as plt9
# import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# matplotlib.use('Agg')
# matplotlib.style.use('ggplot')


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def plot_unity(xdata, ydata, **kwargs):
    mn = min(xdata.min(), ydata.min())
    mx = max(xdata.max(), ydata.max())
    points = np.linspace(mn, mx, 100)
    plt.gca().plot(
        points,
        points,
        color='k',
        marker=None,
        linestyle='--',
        linewidth=1.0
    )


def plot_qq():
    """
    Inspired by https://www.cureffi.org/2012/08/15/qq-plots-with-matplotlib/
    """
    pass


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Compare results from differential expression tests.
            """
    )

    parser.add_argument(
        '-df', '--dataframe',
        action='store',
        dest='df',
        required=True,
        help='Results dataframe. Required columns: \
            gene, pval, log2fc, cell_label_analysed.'
    )

    parser.add_argument(
        '-compare_columns', '--columns_to_compare',
        action='store',
        dest='columns_to_compare',
        default='de_method,covariates',
        help='Comma separated list of columns to use to aggregate results to \
            compare.'
    )

    parser.add_argument(
        '-expr_filt', '--mean_expression_filter',
        action='store',
        dest='mean_expression_filter',
        default=0.1,
        type=float,
        help='Filter genes where mean expression is <= this value. \
            For lowly expressed genes, outliers tend to pull the regression \
            fit and result in false positives. \
            If 0.0 then no filtering is performed.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='output_file',
        default='de_results_comparison',
        help='Basename of output file.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get the parameters
    file = options.df
    cols_compare = options.columns_to_compare.split(',')
    output_file = options.output_file
    mean_expression_filter = options.mean_expression_filter

    # For development
    # cols_compare = [
    #     'condition',
    #     'de_method',
    #     'covariates'
    # ]

    df = pd.read_csv(file, sep='\t')

    # Filter out tests based on mean expression
    if mean_expression_filter != 0.0:
        df = df.loc[(df['mean'] > mean_expression_filter), :]
        output_file = '{}-expression_filter__{}'.format(
            output_file,
            str(mean_expression_filter).replace('.', 'pt')
        )

    # Add a key based on columns to compare (e.g., differential expression
    # method, condition, covariates)
    facet_column = 'facet_key'
    df[facet_column] = df[cols_compare].apply(
        lambda row: '\n'.join(row.values.astype(str)),
        axis=1
    )

    # Add neglog10 value
    small_value = np.empty(len(df), dtype=np.float)
    small_value.fill(0)
    filt = df['pval'] == 0.0
    small_value[filt] = np.nanmin(df['pval'][np.invert(filt)])  # ** 1.5
    df['pval_neglog10'] = np.log10(df['pval'] + small_value) * -1
    df['pval_signedneglog10'] = df['pval_neglog10'] * np.sign(df['log2fc'])

    # For each combination of columns...
    # 1. Plot p-value distribution
    # 2. Plot qq plot
    # 3. Plot -log10 pvalue
    # 4. Plot signed -log10 pvalue

    # Plot p-value distribution
    gplt_pvd = plt9.ggplot(df, plt9.aes(
        x='pval'
    ))
    gplt_pvd = gplt_pvd + plt9.geom_histogram()
    gplt_pvd = gplt_pvd + plt9.theme_bw()
    gplt_pvd = gplt_pvd + plt9.labs(
        x='p-value',
        title=''
    )
    gplt = gplt_pvd + plt9.facet_wrap(
        '~ {}'.format(facet_column),
        ncol=5
    )
    gplt = gplt + plt9.theme(
        strip_text=plt9.element_text(size=5),
        axis_text_x=plt9.element_text(angle=-45, hjust=0)
    )
    n_facets = df[facet_column].nunique()
    gplt.save(
        '{}-pvalue_dist.png'.format(output_file),
        dpi=300,
        width=6*(n_facets/4),
        height=4*(n_facets/4),
        limitsize=False
    )
    # Facet by celltype label
    gplt = gplt_pvd + plt9.facet_grid(
        '{} ~ {}'.format('cell_label_analysed', facet_column)
    )
    gplt = gplt + plt9.theme(
        strip_text=plt9.element_text(size=5),
        axis_text_x=plt9.element_text(angle=-45, hjust=0)
    )
    n_facets_x = df[facet_column].nunique()
    n_facets_y = df['cell_label_analysed'].nunique()
    gplt.save(
        '{}-pvalue_dist-split_cell_label.png'.format(output_file),
        dpi=300,
        width=6*(n_facets_x/4),
        height=4*(n_facets_y/4),
        limitsize=False
    )

    # TODO: QQ plot... color by facet_column with no facet, then add facet
    # for cell_label_analysed

    # Prep to make dataframe to compare p-values across models
    df['row_index'] = df[['gene', 'cell_label_analysed']].apply(
        lambda row: '::'.join(row.values.astype(str)),
        axis=1
    )

    # Make a dataframe to compare p-values across models
    df_plt = df.pivot(
        index='row_index',
        columns='facet_key',
        values='pval_neglog10'
    )
    # Compare -log10 p-value
    sns_plt = sns.pairplot(
        df_plt,
        markers='+',
        kind='reg',
        plot_kws={
            'line_kws': {
                'color': '#ffb90f'
            },
            'scatter_kws': {
                'alpha': 0.1,
                'edgecolor': 'b',
                'linewidth': 1.0
            }
        }
    )
    plt.title('-log10(p-value)')
    sns_plt.map_offdiag(plot_unity)
    sns_plt.savefig(
        '{}-pvalue_neglog10.png'.format(output_file)
    )
    plt.close('all')
    df_plt = df_plt.reset_index()
    df_plt[['gene', 'cell_label_analysed']] = df_plt['row_index'].str.split(
        '::',
        expand=True
    )
    del df_plt['gene']
    df_plt = df_plt.groupby('cell_label_analysed')
    for group_name, df_group in df_plt:
        del df_group['cell_label_analysed']
        sns_plt = sns.pairplot(
            df_group,
            markers='+',
            kind='reg',
            plot_kws={
                'line_kws': {
                    'color': '#ffb90f'
                },
                'scatter_kws': {
                    'alpha': 0.1,
                    'edgecolor': 'b',
                    'linewidth': 1.0
                }
            }
        )
        plt.title('-log10(p-value); cell label: {}'.format(group_name))
        sns_plt.map_offdiag(plot_unity)
        sns_plt.savefig(
            '{}-pvalue_neglog10-cell_label__{}.png'.format(
                output_file,
                group_name
            )
        )
        plt.close('all')

    # Make a dataframe to compare signed p-values across models
    df_plt = df.pivot(
        index='row_index',
        columns='facet_key',
        values='pval_signedneglog10'
    )
    # Compare -log10 p-value
    sns_plt = sns.pairplot(
        df_plt,
        markers='+',
        kind='reg',
        plot_kws={
            'line_kws': {
                'color': '#ffb90f'
            },
            'scatter_kws': {
                'alpha': 0.1,
                'edgecolor': 'b',
                'linewidth': 1.0
            }
        }
    )
    plt.title('Signed -log10(p-value)')
    sns_plt.map_offdiag(plot_unity)
    sns_plt.savefig(
        '{}-pvalue_signedneglog10.png'.format(output_file)
    )
    plt.close('all')
    df_plt = df_plt.reset_index()
    df_plt[['gene', 'cell_label_analysed']] = df_plt['row_index'].str.split(
        '::',
        expand=True
    )
    del df_plt['gene']
    df_plt = df_plt.groupby('cell_label_analysed')
    for group_name, df_group in df_plt:
        del df_group['cell_label_analysed']
        sns_plt = sns.pairplot(
            df_group,
            markers='+',
            kind='reg',
            plot_kws={
                'line_kws': {
                    'color': '#ffb90f'
                },
                'scatter_kws': {
                    'alpha': 0.1,
                    'edgecolor': 'b',
                    'linewidth': 1.0
                }
            }
        )
        plt.title('Signed -log10(p-value); cell label: {}'.format(group_name))
        sns_plt.map_offdiag(plot_unity)
        sns_plt.savefig(
            '{}-pvalue_signedneglog10-cell_label__{}.png'.format(
                output_file,
                group_name
            )
        )
        plt.close('all')


if __name__ == '__main__':
    main()
