#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import scanpy as sc
import plotnine as plt9
# altair: good python plotting package
# import altair
# Seaborn imports below - another good python plotting package
# import seaborn as sns
# import matplotlib
# import matplotlib.pyplot as plt


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def plot_umi_ngene_mt(
    df_plot,
    output_file='plot_umi_ngene_mt',
    facet_column='none',
    color_var='pct_counts_mito_gene'
):
    """Plot plot_umi_ngene_mt to png.

    Parameters
    ----------
    df_plot : pandas.DataFrame
        DataFrame with the followig keys 'total_counts', 'n_genes_by_counts',
        'pct_counts_mito_gene'.
    output_file : string
        Basename of output file.
    facet_column : string
        Column to facet the output by.

    Returns
    -------
    NULL
    """
    if color_var == 'pct_counts_mito_gene':
        color_title = '% MT\n'
    elif color_var == 'cell_passes_qc':
        color_title = 'Cell passed QC\n'
    elif color_var == 'density':
        color_title = 'Density\n'
    else:
        color_title = color_var
    gplt = plt9.ggplot(df_plot, plt9.aes(
        x='total_counts',
        y='n_genes_by_counts',
        color=color_var
    ))
    gplt = gplt + plt9.theme_bw()
    if color_var == 'density':
        gplt = gplt + plt9.geom_density_2d_filled(alpha=0.5)
    gplt = gplt + plt9.geom_point(alpha=0.5, size=0.8)
    gplt = gplt + plt9.scale_x_continuous(
        trans='log10',
        labels=comma_labels,
        minor_breaks=0
    )
    gplt = gplt + plt9.scale_y_continuous(
        trans='log10',
        labels=comma_labels,
        minor_breaks=0
    )
    if color_var == 'pct_counts_mito_gene':
        gplt = gplt + plt9.scale_color_gradient2(
            low='#3B9AB2',
            mid='#EBCC2A',
            high='#F21A00',
            midpoint=50,
            limits=[0, 100]
        )
    elif color_var == 'cell_passes_qc':
        gplt = gplt + plt9.scale_colour_brewer(type='qual', palette='Dark2')
    gplt = gplt + plt9.labs(
        x='Number of molecules',
        y='Number of genes detected',
        color=color_title,
        title=''
    )
    gplt = gplt + plt9.guides(color=plt9.guide_colorbar(ticks=False))
    if facet_column != 'none':
        gplt = gplt + plt9.facet_wrap(
            '~ {}'.format(facet_column),
            ncol=5
        )
        n_samples = df_plot[facet_column].nunique()
        gplt.save(
            '{}.png'.format(output_file),
            dpi=300,
            width=4*(n_samples/4),
            height=4*(n_samples/4),
            limitsize=False
        )
    else:
        gplt.save(
            '{}.png'.format(output_file),
            dpi=300,
            width=4,
            height=4
        )

    # Same plot as above but using seaborn.
    # NOTE: I was having issues setting color scale break points and eventually
    #       gave up

    # Set basic seaborn config.
    # rc = {
    #     'font.size': 12, 'axes.labelsize': 12, 'legend.fontsize': 12,
    #     'axes.titlesize': 12, 'xtick.labelsize': 10, 'ytick.labelsize': 10
    # }
    # sns.set(rc=rc, style='whitegrid')
    #
    # # Set color palette.
    # palette_zissou1 = ["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"]
    # sns.set_palette(sns.color_palette(palette_zissou1))
    #
    # fig, ax = plt.subplots(figsize=(3, 3.5))
    # df_plt['% MT'] = df_plt['pct_counts_mito_gene']
    # ax = sns.scatterplot(
    #     data=df_plt,
    #     ax=ax,
    #     x="total_counts",
    #     y="n_genes_by_counts",
    #     s=5,  # point size
    #     alpha=.5,
    #     hue="% MT",
    #     palette=sns.cubehelix_palette(light=.8, n_colors=6),
    #     linewidth=0
    # )
    # ax.set(xlabel='Number of molecules', ylabel='Number of genes detected')
    # ax.set_xscale('log', basex=10)
    # ax.set_yscale('log', basey=10)
    # for axis in [ax.xaxis, ax.yaxis]:
    #     # Option 1:
    #     # axis.set_major_formatter(
    #     #     matplotlib.ticker.FuncFormatter(
    #     #         lambda y, pos: ('{{:.{:1d}f}}'.format(
    #     #             int(np.maximum(-np.log10(y), 0))
    #     #         )).format(y)
    #     #     )
    #     # )
    #     # Option 2:
    #     # axis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #     # Option 3:
    #     axis.set_major_formatter(
    #         matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    #     )
    # ax.legend(
    #     frameon=False,
    #     loc='center left',
    #     bbox_to_anchor=(1, 0.5),
    #     ncol=1,
    #     fontsize=10,
    # )
    # fig.savefig(
    #     '{}.png'.format('test'),
    #     dpi=300,
    #     bbox_inches='tight'
    # )
    # plt.close(fig)  # Close the figure.
    #
    # # Seaborn fact plot.
    # g = sns.FacetGrid(df_plot, col="sanger_sample_id")
    # g = g.map(plt.hist, "n_genes_by_counts")
    # g.savefig(
    #     '{}.png'.format('test'),
    #     dpi=300,
    #     bbox_inches='tight'
    # )
    # plt.close(g)   # Close the figure.


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter and merge 10x data. Save to AnnData object.
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
        '-of', '--output_file',
        action='store',
        dest='of',
        default='plot_umi_ngene_mt',
        help='Basename of output png file. Will have .png appended.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--facet_columns',
        action='store',
        dest='facet_columns',
        default='none',
        help='Column to facet plot by.\
            (default: No facet.)'
    )

    options = parser.parse_args()

    # Scanpy settings
    # sc.settings.figdir = os.getcwd()  # figure output directory 2 match base.
    # # sc.settings.n_jobs = options.ncpu  # number CPUs
    # # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save=300)

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)

    # Get a list of the facets to plot.
    facets_to_plot = options.facet_columns.split(',')
    if len(facets_to_plot) == 0:
        facets_to_plot = ['none']

    # Plot the data.
    for facet in facets_to_plot:
        plot_umi_ngene_mt(
            df_plot=adata.obs,
            output_file='plot_umi_ngene_mt.facet={}-{}'.format(
                facet,
                options.of
            ),
            facet_column=facet,
            color_var='pct_counts_mito_gene'
        )
        plot_umi_ngene_mt(
            df_plot=adata.obs,
            output_file='plot_umi_ngene_density.facet={}-{}'.format(
                facet,
                options.of
            ),
            facet_column=facet,
            color_var='density'
        )
        if 'cell_passes_qc' in adata.obs:
            plot_umi_ngene_mt(
                df_plot=adata.obs,
                output_file='plot_umi_ngene_cellpassqc.facet={}-{}'.format(
                    facet,
                    options.of
                ),
                facet_column=facet,
                color_var='cell_passes_qc'
            )


if __name__ == '__main__':
    main()
