#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import csv
import numpy as np
import pandas as pd
import scrublet as scr
import skimage.filters as skif
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import plotnine as plt9
# from statannot import add_stat_annotation

sns.set(style='whitegrid')


# TODO: make better UMAP
# TODO: make better plot of z scores


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def plot_scrub_hist(
    scrub,
    threshold,
    output_file='scrublet_histogram',
    scale_y_log10=False,
    density=False,
    zscores=False
):
    # Make a pandas dataframe of the data
    if not zscores:
        df_tmp1 = pd.DataFrame(
            data=scrub.doublet_scores_obs_,
            columns=['scores']
        )
        df_tmp2 = pd.DataFrame(
            data=scrub.doublet_scores_sim_,
            columns=['scores']
        )
        x_axis_label = 'Doublet score'
    else:
        # Calculate z scores as in the scrublet source code (see call_doublets
        # function).
        df_tmp1 = pd.DataFrame(
            data=scrub.z_scores_,
            columns=['scores']
        )
        zscore_sim = (
            scrub.doublet_scores_sim_ - threshold
        ) / scrub.doublet_errors_sim_
        df_tmp2 = pd.DataFrame(
            data=zscore_sim,
            columns=['scores']
        )
        x_axis_label = 'Doublet zscore'

    df_tmp1['type'] = 'Observed'
    df_tmp2['type'] = 'Simulated'

    df_plt = pd.concat([df_tmp1, df_tmp2])

    # Plot the data usig plotnine
    gplt = plt9.ggplot(df_plt, plt9.aes(
        x='scores'
    ))
    gplt = gplt + plt9.theme_bw()
    if density:
        gplt = gplt + plt9.geom_density(alpha=0.8)
    else:
        gplt = gplt + plt9.geom_histogram(alpha=0.8)
    if not zscores:
        gplt = gplt + plt9.geom_vline(xintercept=threshold, linetype='dotted')
    if scale_y_log10:
        gplt = gplt + plt9.scale_y_continuous(
            trans='log10',
            #labels=comma_labels,
            minor_breaks=0
        )
    gplt = gplt + plt9.labs(
        x=x_axis_label,
        # y='Counts',
        title=''
    )
    gplt = gplt + plt9.facet_wrap(
        '~ {}'.format('type'),
        scales='free_y'
    )

    # Sort out whitespace issue in plotnine:
    # https://github.com/has2k1/plotnine/issues/106
    gplt = gplt + plt9.theme(subplots_adjust={'wspace': 0.35})

    gplt.save(
        '{}.png'.format(output_file),
        dpi=300,
        width=4,
        height=2,
        limitsize=False
    )

# def scanpy_merge(
#     tenx_data,
#     metadata,
#     metadata_key,
#     output_file,
#     params_dict=dict()
# ):
#     """Merge 10x data.
#
#     Parameters
#     ----------
#     tenx_data : pandas.DataFrame
#         Description of parameter `tenx_data`.
#     metadata : pandas.DataFrame
#         Description of parameter `metadata`.
#     output_file : string
#         Description of parameter `output_file`.
#     metadata_key : string
#         Column in metadata that matches the "experiment_id" column in
#         tenx_data (the default is "sanger_sample_id").
#
#     Returns
#     -------
#     output_file : string
#         output_file
#     """
#     # check the tenx_data
#     # check for required columns
#     tenx_data_required_cols = set(['experiment_id', 'data_path_10x_format'])
#     if not tenx_data_required_cols.issubset(tenx_data.columns):
#         raise Exception('Invalid tenx_data.')
#     # check no duplicate sample ids
#     vals, counts = np.unique(tenx_data['experiment_id'], return_counts=True)
#     if np.sum(counts > 1):
#         raise Exception('Error {} duplicate experiment_ids:\t{}'.format(
#             np.sum(counts > 1),
#             np.array2string(vals[counts > 1])
#         ))
#     # check all files exist
#     filt = tenx_data['data_path_10x_format'].apply(
#         lambda x: os.path.exists(x)
#     )
#     if np.sum(filt > 1):
#         raise Exception('Error {} data_path_10x_format missing:\t{}'.format(
#             np.sum(filt > 1),
#             np.array2string(tenx_data['data_path_10x_format'][filt])
#         ))
#
#     # Init default values for params_dict
#     params_filters_check = [
#         'cell_filters',
#         'downsample_cells_fraction',
#         'downsample_cells_n',
#         'downsample_feature_counts'
#     ]
#     for i in params_filters_check:
#         if i not in params_dict:
#             if i == 'cell_filters':
#                 params_dict[i] = {'value': []}
#             else:
#                 params_dict[i] = {'value': ''}
#     # Check for validity of filters specified in params_dict
#     param_filters_check = [
#         'downsample_cells_fraction',
#         'downsample_cells_n'
#     ]
#     if all(params_dict[k]['value'] != '' for k in param_filters_check):
#         raise Exception(
#             'Error check the params. Both {} and {} are set.'.format(
#                 'downsample_cells_fraction',
#                 'downsample_cells_n'
#             )
#         )
#
#     # iterate over samples and load data
#     adatasets = []
#     n_adatasets = 1
#     for idx, row in tenx_data.iterrows():
#         # load the data
#         adata = sc.read_10x_mtx(
#             path=row['data_path_10x_format'],
#             # var_names='gene_symbols',
#             var_names='gene_ids',
#             make_unique=False
#         )
#
#         adata = check_adata(adata, row['experiment_id'])
#
#         # label mitochondrial genes
#         # mito_gene_list = sc.queries.mitochondrial_genes() # another query
#         adata.var['mito_gene'] = [
#             x.startswith('MT-') for x in adata.var['gene_symbols']
#         ]
#         # use this if var_names='gene_symbols' in sc.read_10x_mtx
#         # adata.var['mito_gene'] = [
#         #     x.startswith('MT-') for x in adata.var_names
#         # ]
#
#         # Add in sample metadata.
#         # NOTE: it would be more memory efficient to stash this in
#         #       unstructured dict-like annotation (adata.uns)
#         metadata_smpl = metadata[
#             metadata[metadata_key] == row['experiment_id']
#         ]
#         for col in metadata_smpl.columns:
#             adata.obs[col] = np.repeat(metadata_smpl[col].values, adata.n_obs)
#
#         # Calculate basic qc metrics for this sample.
#         # NOTE: n_genes_by_counts == number of genes with > 0 counts
#         #       adata.obs['n_genes'] = (adata.X > 0).sum(axis = 1) is same as
#         #       adata.obs['n_genes_by_counts']
#         vars_prior_metrics = adata.var_keys()
#         sc.pp.calculate_qc_metrics(adata, qc_vars=['mito_gene'], inplace=True)
#
#         # Apply cell filter.
#         # adata = adata[selected_cells, :]
#         # Apply gene filter
#         # adata = adata[:, selected_genes]
#
#         # Apply downsampling if needed.
#         if params_dict['downsample_cells_fraction']['value'] != '':
#             n_cells_start = adata.n_obs
#             sc.pp.subsample(
#                 adata,
#                 fraction=float(
#                     params_dict['downsample_cells_fraction']['value']
#                 ),
#                 copy=False
#             )
#             print('[{}] cell downsample applied: {} dropped {} remain'.format(
#                 row['experiment_id'],
#                 n_cells_start - adata.n_obs,
#                 adata.n_obs
#             ))
#         elif params_dict['downsample_cells_n']['value'] != '':
#             n_cells_start = adata.n_obs
#             sc.pp.subsample(
#                 adata,
#                 n_obs=int(params_dict['downsample_cells_n']['value']),
#                 copy=False
#             )
#             print('[{}] cell downsample applied: {} dropped {} remain'.format(
#                 row['experiment_id'],
#                 n_cells_start - adata.n_obs,
#                 adata.n_obs
#             ))
#         if params_dict['downsample_feature_counts']['value'] != '':
#             fraction = params_dict['downsample_feature_counts']['value']
#             target_counts_per_cell = adata.obs['total_counts'].apply(
#                 lambda x: int(x * fraction)
#             ).values
#             sc.pp.downsample_counts(
#                 adata,
#                 counts_per_cell=target_counts_per_cell
#             )
#         # Apply cell QC filters.
#         if len(params_dict['cell_filters']['value']) > 0:
#             n_cells_start = adata.n_obs
#             for filter_query in params_dict['cell_filters']['value']:
#                 if filter_query != '':
#                     adata = adata[adata.obs.query(filter_query).index, :]
#                     print('[{}] {} "{}": {} dropped {} remain'.format(
#                         'cell QC applied',
#                         row['experiment_id'],
#                         filter_query,
#                         n_cells_start - adata.n_obs,
#                         adata.n_obs
#                     ))
#             print('[{}] after all cell QC: {} dropped {} remain'.format(
#                 row['experiment_id'],
#                 n_cells_start - adata.n_obs,
#                 adata.n_obs
#             ))
#
#         # Print the number of cells and genes for this sample.
#         print('[{}] {} obs (cells), {} vars (genes)'.format(
#             row['experiment_id'],
#             adata.n_obs,
#             adata.n_vars
#         ))
#
#         # Comment code below to keep the vars (gene) output from
#         # calculate_qc_metrics *per sample*. If we do this, then in
#         # adata_merged.var, we will have duplicated # measures according to
#         # each sample (e.g., n_cells_by_counts-0, # n_cells_by_counts-1,
#         # n_cells_by_counts-3).
#         #
#         # Code below removes such output.
#         adata.var = adata.var[vars_prior_metrics]
#
#         # If we still have cells after filters, add to our list of data.
#         if adata.n_obs > 0:
#             adatasets.append(adata)
#             n_adatasets += 1
#         else:
#             raise Exception(
#                 'Error invalid tenx_data file. Missing coluns.'
#             )
#
#     # Merge all of the data together.
#     adata_merged = adatasets[0].concatenate(*adatasets[1:])
#     adata_merged = check_adata(adata_merged, 'adata_merged')
#
#     # Re-calculate basic qc metrics for the whole dataset.
#     obs_prior = adata_merged.obs.copy()
#     sc.pp.calculate_qc_metrics(
#         adata_merged,
#         qc_vars=['mito_gene'],
#         inplace=True
#     )
#     adata_merged.obs = obs_prior
#
#     # Possible additional basic filtering on the full dataset.
#     # sc.pp.filter_cells(adata, min_genes=200)
#     # sc.pp.filter_genes(adata, min_cells=1)
#
#     print('[adata_merged] {} obs, {} vars'.format(
#         adata_merged.n_obs,
#         adata_merged.n_vars
#     ))
#
#     # output_file = output_dir + "/adata"
#     adata_merged.write('{}.h5ad'.format(output_file), compression='gzip')
#     # adata_merged.write_csvs(output_file)
#     # adata_merged.write_loom(output_file+".loom")
#
#     return(output_file)


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
        '-txd', '--tenxdata_dir',
        action='store',
        dest='txd',
        default='',
        required=False,
        help='Path to directory with data in 10x matrix format.'
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
        '-emr', '--expected_multiplet_rate',
        action='store',
        dest='emr',
        default=0.0,
        type=float,
        help='Expected multiplet rate. If 0 then automatically predict doublet\
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
    verbose = True

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
    if options.emr == 0.0:
        multiplet_rate = 0.0007589 * cells_recovered + 0.0527214
        multiplet_rate = multiplet_rate / 100.0
    else:
        multiplet_rate = options.emr
    # expected_n_doublets = multiplet_rate * cells_recovered
    if verbose:
        print('multiplet_rate:\t{}'.format(multiplet_rate))

    n_simulated_doublets = 100000.0
    #n_simulated_doublets = 20000.0
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
        sim_doublet_ratio=n_simulated_doublets/len(adata),  # Default is 2.0
        expected_doublet_rate=multiplet_rate  # Default is 0.1
    )

    # Run the scrublet pipeline with default parameters.
    # This function performs:
    # * Doublet simulation
    # * Normalization, gene filtering, rescaling, PCA
    # * Doublet score calculation
    # * Doublet score threshold detection and doublet calling
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        verbose=False
    )

    # Calculate the threshold for calling doublets
    # The default method for scrublet is `threshold_minimum`, but in our
    # hands it looks like `threshold_otsu` looks better.
    # threshold = skif.threshold_minimum(scrub.doublet_scores_sim_)
    #
    # Estimate the number of bins using Freedman–Diaconis rule which is robust
    # to large datasets.
    hist, _ = np.histogram(scrub.doublet_scores_sim_, bins='fd')
    n_bins = len(hist)
    threshold = skif.threshold_otsu(scrub.doublet_scores_sim_, nbins=n_bins)
    print('threshold_minimum:\t{}'.format(
        skif.threshold_minimum(scrub.doublet_scores_sim_, nbins=n_bins)
    ))
    print('threshold_triangle:\t{}'.format(
        skif.threshold_triangle(scrub.doublet_scores_sim_, nbins=n_bins)
    ))
    print('threshold_yen:\t{}'.format(
        skif.threshold_yen(scrub.doublet_scores_sim_, nbins=n_bins)
    ))
    print('threshold_otsu:\t{}'.format(threshold))
    # Call doublets using the otsu threshold
    predicted_doublets = scrub.call_doublets(
        threshold=threshold,
        verbose=verbose
    )
    print(
        'If automatic threshold is poor, adjust threshold with',
        'scrub.call_doublets(threshold=<my_custom_threshold>)'
    )
    adata.obs['scrublet__predicted_doublets'] = scrub.predicted_doublets_
    adata.obs['scrublet__doublet_scores'] = scrub.doublet_scores_obs_
    adata.obs['scrublet__doublet_zscores'] = scrub.z_scores_

    # Get the number of doublets one would expect if the 10x prediction were
    # spot on. Taken from
    # https://github.com/vib-singlecell-nf/scrublet/blob/master/bin/sc_doublet_detection.py
    #
    # Estimate the multiplet rate using coefficients from a lm predicting
    # multiplet rate from recovered cells in data distributed by 10x for
    # their 3.1 chemistry.
    # https://support.10xgenomics.com/single-cell-gene-expression/library-prep/doc/user-guide-chromium-single-cell-3-reagent-kits-user-guide-v31-chemistry
    # cells_recovered = len(adata)
    # multiplet_rate = 0.0007589 * cells_recovered + 0.0527214
    # multiplet_rate = multiplet_rate / 100.0
    # doublet_cells = adata.obs['scrublet__doublet_scores'].sort_values(
    #     ascending=False
    # ).head(
    #     n=expected_n_doublets
    # ).index
    # adata.obs['scrublet__predicted_doublets_based_10x_doublet_rate'] = False
    # adata.obs.loc[
    #     doublet_cells,
    #     'scrublet__predicted_doublets_based_10x_doublet_rate'
    # ] = True

    # Save the results.
    cols_save = [
        'scrublet__doublet_scores',
        'scrublet__predicted_doublets',
        'scrublet__doublet_zscores'
    ]
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
    plot_scrub_hist(
        scrub=scrub,
        threshold=threshold,
        scale_y_log10=False,
        output_file='{}-histogram_doublet_scores_log.pdf'.format(out_file_base)
    )
    plot_scrub_hist(
        scrub=scrub,
        threshold=threshold,
        scale_y_log10=True,
        output_file='{}-histogram_doublet_scores.pdf'.format(out_file_base)
    )
    plot_scrub_hist(
        scrub=scrub,
        threshold=threshold,
        density=True,
        scale_y_log10=False,
        output_file='{}-density_doublet_scores.pdf'.format(out_file_base)
    )
    plot_scrub_hist(
        scrub=scrub,
        threshold=threshold,
        zscores=True,
        output_file='{}-histogram_doublet_zscores.pdf'.format(out_file_base)
    )
    fig, ax = scrub.plot_histogram(
        scale_hist_obs='linear',
        scale_hist_sim='linear'
    )
    fig.savefig(
        '{}-histogram_doublet_scores_v2.pdf'.format(out_file_base),
        dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)

    # Plot UMAP embedding.
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_))
    fig, ax = scrub.plot_embedding('UMAP', order_points=True)
    fig.savefig(
        '{}-umap_doublet_scores.png'.format(out_file_base),
        dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)

    fig, ax = scrub.plot_embedding('UMAP', score='zscore', order_points=True)
    fig.savefig(
        '{}-umap_doublet_zscores.png'.format(out_file_base),
        dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)

    # Plot tSNE embedding.
    # scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))

    # Plot the average number of UMIs in multiplets vs singlets.
    if 'total_counts' not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
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
        '{}-boxplot_total_umi_counts.png'.format(out_file_base),
        dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)  # Close the figure.


if __name__ == '__main__':
    main()
