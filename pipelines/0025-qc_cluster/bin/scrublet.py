#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import scanpy as sc
import csv
import scrublet as scr
import matplotlib.pyplot as plt
import seaborn as sns
# from statannot import add_stat_annotation

sns.set(style='whitegrid')


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read 10x data or AnnData object and annotates scrublets.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-txd', '--tenxdata_file',
        action='store',
        dest='txd',
        default='',
        required=False,
        help='File with the following headers: experiment_id\
            data_path_10x_format.'
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        default='',
        required=False,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-edr', '--expected_doublet_rate',
        action='store',
        dest='edr',
        default=0,
        type=float,
        help='Expected doublet rate. If 0 then automatically predict doublet\
            rate based on the number of recovered cells using published\
            doublet rates from 10x v3.1 chemistry.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='scrublet',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Fixed settings.
    # verbose = True

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    # sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    # Get the out file base.
    out_file_base = options.of

    # Load the data.
    if options.txd != '':
        adata = sc.read_10x_mtx(
            path=options.txd,
            # var_names='gene_symbols',
            var_names='gene_ids',
            make_unique=False
        )
    elif options.h5 != '':
        adata = sc.read_h5ad(filename=options.h5)
        print(
            'Scrublet uses counts. Assuming adata.X are counts'
        )
    else:
        raise Exception(
            'Error invalid input.'
        )

    # If no expected doublet rate, then estimate the doublet rate using
    # the coefficients from a lm predicting doublet rate from recovered cells
    # in data distributed by 10x for their 3.1 chemistry.
    # https://support.10xgenomics.com/single-cell-gene-expression/library-prep/doc/user-guide-chromium-single-cell-3-reagent-kits-user-guide-v31-chemistry
    cells_recovered = len(adata)
    if options.edr == 0.0:
        multiplet_rate = 0.0007589 * cells_recovered + 0.0527214
    else:
        multiplet_rate = options.edr
    expected_n_doublets = multiplet_rate / 100.0 * cells_recovered

    # How to set expected_doublet_rate.
    # The authors note that the method is not that sensitive to this parameter
    # https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb
    # From Chromium Single Cell 3' Reagent Kits User Guide (v2 Chemistry):
    #   https://support.10xgenomics.com/permalink/3vzDu3zQjY0o2AqkkkI4CC
    #   the expected doublet rate is 0.06.
    # From Chromium Single Cell 3' Reagent Kits User Guide (v3.1 Chemistry):
    #   https://support.10xgenomics.com/single-cell-gene-expression/library-prep/doc/user-guide-chromium-single-cell-3-reagent-kits-user-guide-v31-chemistry
    #   the expected doublet rate is ~3.9% for ~8000 input cells and ~1.6% for
    #   ~3,200 input cells.
    scrub = scr.Scrublet(
        counts_matrix=adata.X,
        sim_doublet_ratio=2.0,  # Default is 2.0
        expected_doublet_rate=multiplet_rate  # Default is 0.1
    )

    # Run the scrublet pipeline with default parameters.
    # This function performs:
    # * Doublet simulation
    # * Normalization, gene filtering, rescaling, PCA
    # * Doublet score calculation
    # * Doublet score threshold detection and doublet calling
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        verbose=True
    )
    adata.obs['scrublet__doublet_scores'] = doublet_scores
    adata.obs['scrublet__predicted_doublets'] = predicted_doublets

    # Get the number of doublets one would expect if the 10x prediction were
    # spot on.
    doublet_cells = adata.obs['scrublet__doublet_scores'].sort_values(
        ascending=False
    ).head(
        n=expected_n_doublets
    ).index
    adata.obs['scrublet__predicted_doublets_based_10x_doublet_rate'] = False
    adata.obs.loc[
        doublet_cells,
        'scrublet__predicted_doublets_based_10x_doublet_rate'
    ] = True

    # Save the results.
    cols_save = ['scrublet__doublet_scores', 'scrublet__predicted_doublets']
    adata.obs[cols_save].to_csv(
        '{}-scrublet.tsv.gz'.format(out_file_base),
        sep='\t',
        index=True,
        quoting=csv.QUOTE_NONNUMERIC,
        index_label='cell_barcode',
        na_rep='',
        compression='gzip'
    )

    # Plot a histogram of doublets.
    fig, ax = scrub.plot_histogram(
        scale_hist_obs='linear',
        scale_hist_sim='linear'
    )
    fig.savefig(
        '{}-histogram_doublet_scores.pdf'.format(out_file_base),
        dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)
    fig, ax = scrub.plot_histogram(
        scale_hist_obs='log',
        scale_hist_sim='log'
    )
    fig.savefig(
        '{}-histogram_doublet_scores_log.pdf'.format(out_file_base),
        dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)

    print(
        'If automatic threshold is poor, adjust threshold with',
        'scrub.call_doublets(threshold=<my_custom_threshold>)'
    )

    # How to call doublets manually if the default threshold does not look
    # great.
    # scrub.call_doublets(threshold=0.25)

    # Plot the average number of UMIs in multiplets vs singlets.
    fig, ax = plt.subplots(figsize=(3, 3.5))
    ax = sns.boxplot(
        data=adata.obs,
        x='scrublet__predicted_doublets',
        y='total_counts'
        # hue='scrublet__predicted_doublets'
    )
    ax.set_yscale('log', basey=10)
    ax.set(xlabel='Predicted multiplet', ylabel='Number of molecules')
    # NOTE: I could not get p-value annotation to work.
    # ax, test_results = add_stat_annotation(
    #     ax,
    #     data=adata.obs,
    #     x='scrublet__predicted_doublets',
    #     y='total_counts',
    #     #hue='scrublet__predicted_doublets',
    #     box_pairs=[('True', 'False')],
    #     test='Mann-Whitney',
    #     text_format='full',
    #     loc='inside'
    # )
    # plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))
    fig.savefig(
        '{}.png'.format('test'),
        dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)  # Close the figure.

    # Plot UMAP embedding.
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_))
    fig, ax = scrub.plot_embedding('UMAP', order_points=True)
    fig.savefig(
        '{}-umap_doublet_scores_log.png'.format(out_file_base),
        dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)

    # Plot tSNE embedding.
    # scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))


if __name__ == '__main__':
    main()
