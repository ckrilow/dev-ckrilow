#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import yaml
import numpy as np
import pandas as pd
import scanpy as sc

# sc verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.settings.verbosity = 3
# sc.logging.print_versions()
# sc.settings.set_figure_params(dpi=80)


def check_adata(adata, adata_id):
    """Check adata."""
    # get var names (genes) that are not unique
    vals, counts = np.unique(adata.var_names, return_counts=True)
    if np.sum(counts > 1):
        print('[{}] fixing {} duplicate var_names:\t{}'.format(
            adata_id,
            np.sum(counts > 1),
            np.array2string(vals[counts > 1])
        ))
        adata.var_names_make_unique()

    # get obs_names (cell barcodes) that are not unique
    vals, counts = np.unique(adata.obs_names, return_counts=True)
    if np.sum(counts > 1):
        raise Exception('[{}] error {} duplicate obs_names:\t{}'.format(
            adata_id,
            np.sum(counts > 1),
            np.array2string(vals[counts > 1])
        ))

    return(adata)


def scanpy_merge(
    tenx_data,
    metadata,
    metadata_key,
    output_file,
    params_dict=dict()
):
    """Merge 10x data.

    Parameters
    ----------
    tenx_data : pandas.DataFrame
        Description of parameter `tenx_data`.
    metadata : pandas.DataFrame
        Description of parameter `metadata`.
    output_file : string
        Description of parameter `output_file`.
    metadata_key : string
        Column in metadata that matches the "experiment_id" column in
        tenx_data (the default is "sanger_sample_id").

    Returns
    -------
    output_file : string
        output_file
    """
    # check the tenx_data
    # check for required columns
    tenx_data_required_cols = set(['experiment_id', 'data_path_10x_format'])
    if not tenx_data_required_cols.issubset(tenx_data.columns):
        raise Exception('Invalid tenx_data.')
    # check no duplicate sample ids
    vals, counts = np.unique(tenx_data['experiment_id'], return_counts=True)
    if np.sum(counts > 1):
        raise Exception('Error {} duplicate experiment_ids:\t{}'.format(
            np.sum(counts > 1),
            np.array2string(vals[counts > 1])
        ))
    # check all files exist
    filt = tenx_data['data_path_10x_format'].apply(
        lambda x: os.path.exists(x)
    )
    if np.sum(filt > 1):
        raise Exception('Error {} data_path_10x_format missing:\t{}'.format(
            np.sum(filt > 1),
            np.array2string(tenx_data['data_path_10x_format'][filt])
        ))

    # Init default values for params_dict
    params_filters_check = [
        'cell_filters',
        'downsample_cells_fraction',
        'downsample_cells_n',
        'downsample_feature_counts'
    ]
    for i in params_filters_check:
        if i not in params_dict:
            if i == 'cell_filters':
                params_dict[i] = {'value': []}
            else:
                params_dict[i] = {'value': ''}
    # Check for validity of filters specified in params_dict
    param_filters_check = [
        'downsample_cells_fraction',
        'downsample_cells_n'
    ]
    if all(params_dict[k]['value'] != '' for k in param_filters_check):
        raise Exception(
            'Error check the params. Both {} and {} are set.'.format(
                'downsample_cells_fraction',
                'downsample_cells_n'
            )
        )

    # iterate over samples and load data
    adatasets = []
    n_adatasets = 1
    for idx, row in tenx_data.iterrows():
        # load the data
        adata = sc.read_10x_mtx(
            path=row['data_path_10x_format'],
            # var_names='gene_symbols',
            var_names='gene_ids',
            make_unique=False
        )

        adata = check_adata(adata, row['experiment_id'])

        # label mitochondrial genes
        # mito_gene_list = sc.queries.mitochondrial_genes() # another query
        adata.var['mito_gene'] = [
            x.startswith('MT-') for x in adata.var['gene_symbols']
        ]
        # use this if var_names='gene_symbols' in sc.read_10x_mtx
        # adata.var['mito_gene'] = [
        #     x.startswith('MT-') for x in adata.var_names
        # ]

        # Add in sample metadata.
        # NOTE: it would be more memory efficient to stash this in
        #       unstructured dict-like annotation (adata.uns)
        metadata_smpl = metadata[
            metadata[metadata_key] == row['experiment_id']
        ]
        for col in metadata_smpl.columns:
            adata.obs[col] = np.repeat(metadata_smpl[col].values, adata.n_obs)

        # Calculate basic qc metrics for this sample.
        vars_prior_metrics = adata.var_keys()
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mito_gene'], inplace=True)

        # Apply cell filter.
        # adata = adata[selected_cells, :]
        # apply gene filter
        # adata = adata[:, selected_genes]

        # Apply downsampling if needed.
        if params_dict['downsample_cells_fraction']['value'] != '':
            sc.pp.subsample(
                adata,
                fraction=params_dict['downsample_cells_fraction']['value'],
                copy=False
            )
        elif params_dict['downsample_cells_fraction']['value'] != '':
            sc.pp.subsample(
                adata,
                n_obs=params_dict['downsample_cells_n']['value'],
                copy=False
            )
        if params_dict['downsample_feature_counts']['value'] != '':
            fraction = params_dict['downsample_feature_counts']['value']
            target_counts_per_cell = adata.obs['total_counts'].apply(
                lambda x: int(x * fraction)
            ).values
            sc.pp.downsample_counts(
                adata,
                counts_per_cell=target_counts_per_cell
            )
        # Apply cell QC filters.
        if len(params_dict['cell_filters']['value']) > 0:
            n_cells_start = adata.n_obs
            for filter_query in params_dict['cell_filters']['value']:
                adata = adata[adata.obs.query(filter_query).index, :]
                print('[{}] cell QC applied "{}": {} cells dropped'.format(
                    row['experiment_id'],
                    filter_query,
                    n_cells_start - adata.n_obs
                ))
            print('[{}] cell QC applied: {} total cells dropped'.format(
                row['experiment_id'],
                n_cells_start - adata.n_obs
            ))

        # Print the number of cells and genes for this sample.
        print('[{}] {} obs (cells), {} vars (genes)'.format(
            row['experiment_id'],
            adata.n_obs,
            adata.n_vars
        ))

        # Comment code below to keep the vars (gene) output from
        # calculate_qc_metrics *per sample*. If we do this, then in
        # adata_merged.var, we will have duplicated # measures according to
        # each sample (e.g., n_cells_by_counts-0, # n_cells_by_counts-1,
        # n_cells_by_counts-3).
        #
        # Code below removes such output.
        adata.var = adata.var[vars_prior_metrics]

        # If we still have cells after filters, add to our list of data.
        if adata.n_obs > 0:
            adatasets.append(adata)
            n_adatasets += 1

    # Merge all of the data together.
    adata_merged = adatasets[0].concatenate(*adatasets[1:])
    adata_merged = check_adata(adata_merged, 'adata_merged')

    # Re-calculate basic qc metrics for the whole dataset.
    obs_prior = adata_merged.obs.copy()
    sc.pp.calculate_qc_metrics(
        adata_merged,
        qc_vars=['mito_gene'],
        inplace=True
    )
    adata_merged.obs = obs_prior

    # Possible additional basic filtering on the full dataset.
    # sc.pp.filter_cells(adata, min_genes=200)
    # sc.pp.filter_genes(adata, min_cells=1)

    print('[adata_merged] {} obs, {} vars'.format(
        adata_merged.n_obs,
        adata_merged.n_vars
    ))

    # output_file = output_dir + "/adata"
    adata_merged.write('{}.h5'.format(output_file), compression='gzip')
    # adata_merged.write_csvs(output_file)
    # adata_merged.write_loom(output_file+".loom")

    return(output_file)


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
        '-txd', '--tenxdata_file',
        action='store',
        dest='txd',
        required=True,
        help='File with the following headers: experiment_id\
            data_path_10x_format.'
    )

    parser.add_argument(
        '-mf', '--metadata_file',
        action='store',
        dest='mf',
        required=True,
        help='File with metadata on samples matching sanger_sample_id in\
            tenxdata_file.'
    )

    parser.add_argument(
        '-mcd', '--metadata_columns_delete',
        action='store',
        dest='mcd',
        default='sample_status,study,study_id',
        help='Comma seperated list of columns to delete in metadata_file.\
            If "" then no columns are deleted. Not whitespace should be\
            represented with an underscore (_).\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-mk', '--metadata_key',
        action='store',
        dest='mk',
        default='sanger_sample_id',
        help='Key to link metadata to tenxdata_file experiment_id column.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=2,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-pyml', '--params_yaml',
        action='store',
        dest='pyml',
        default='',
        help='YAML file containing cell filtering and downsampling (e.g.,\
            of total number of cells or reads per cell) parameters.\
            If file is not provided, no filtering or downsampling is\
            performed.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='adata',
        help='Basename of output anndata file, assuming output in current \
            working directory. Will have .h5 appended.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    # NOTE:
    # - Could change yaml params file to include a list per sample and the
    #   filters one wants to use for that sample
    # - Could allow the input to be either 10x dir or AnnData/loom object.
    #   Use AnnData/loom object would be useful if we add doublet scores
    #   or other scores prior to filtering.

    # Read in the paramerters for downsampling and cell filters.
    if options.pyml == '':
        params_dict = {}
    else:
        with open(options.pyml, 'r') as f:
            params_dict = yaml.safe_load(f)
        params_dict = params_dict['sample_qc']

    # Load a file of the samples to analyse
    tenx_data = pd.read_csv(options.txd, sep='\t')
    # tenx_data = tenx_data.rename(columns={
    #     'sanger_sample_id': 'experiment_id',
    #     'file_path': 'data_path_10x_format'
    # })
    tenx_data_check = [
        'experiment_id',
        'data_path_10x_format'
    ]
    if not all(k in tenx_data.columns for k in tenx_data_check):
        raise Exception(
            'Error invalid tenx_data file. Missing coluns.'
        )

    # Load the metadata
    metadata = pd.read_csv(options.mf, sep='\t')
    metadata.columns = metadata.columns.str.strip(
        ).str.replace(' ', '_').str.lower()

    # Delete the metadata columns that we do not want.
    if options.mcd != '':
        for i in options.mcd.split(','):
            if i in metadata.columns:
                metadata = metadata.drop(i, axis=1)

    # Make sure the matching key exists.
    if options.mk not in metadata.columns:
        raise Exception(
            'Error cannot find metadata_key in metadata.'
        )

    # Run the merge function.
    out_file = scanpy_merge(
        tenx_data,
        metadata,
        metadata_key=options.mk,
        output_file=options.of,
        params_dict=params_dict
    )
    print(out_file)


if __name__ == '__main__':
    main()
