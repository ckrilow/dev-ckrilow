#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import csv
import os
import yaml
import random
import warnings
import numpy as np
import pandas as pd
import scanpy as sc

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)

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
    params_dict=dict(),
    cellmetadata_filepaths=None
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
        tenx_data.

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

    # Check the cellmetadata_filepaths if we have it
    if cellmetadata_filepaths is not None:
        # check all files exist
        filt = cellmetadata_filepaths['data_path_cellmetadata'].apply(
            lambda x: os.path.exists(x)
        )
        if np.sum(filt > 1):
            raise Exception(
                'Error {} data_path_cellmetadata missing:\t{}'.format(
                    np.sum(filt > 1),
                    np.array2string(
                        cellmetadata_filepaths['data_path_cellmetadata'][filt]
                    )
                )
            )

    # Init default values for params_dict
    params_filters_check = [
        'cell_filters',
        'downsample_cells_fraction',
        'downsample_cells_n',
        'downsample_feature_counts'
    ]
    for i in params_filters_check:
        if i not in params_dict:
            if i != 'cell_filters':
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

    # Init a dictionary to record the original number of cells per sample and
    # the number of cells after each filter.
    n_cells_dict = {}

    # If true, then filtered cells are dropped prior to merging the data.
    # This will save disk space.
    drop_filtered_cells = False

    # Iterate over samples and load data
    adatasets = []
    n_adatasets = 1
    for idx, row in tenx_data.iterrows():
        # Load the data
        adata = sc.read_10x_mtx(
            path=row['data_path_10x_format'],
            # var_names='gene_symbols',
            var_names='gene_ids',
            make_unique=False
        )

        adata = check_adata(adata, row['experiment_id'])

        # Record the total number of cells for this experiment_id
        n_cells_dict[row['experiment_id']] = {}
        n_cells_dict[row['experiment_id']]['before_filters'] = adata.n_obs

        # Label mitochondrial genes
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

        # Ensure we have experiment_in the final dataframe.
        if 'experiment_id' not in adata.obs.columns:
            adata.obs['experiment_id'] = adata.obs[metadata_key]

        # Add in per cell metadata if we have it.
        if cellmetadata_filepaths is not None:
            if row['experiment_id'] in cellmetadata_filepaths.index:
                cellmetadata = pd.read_csv(
                    cellmetadata_filepaths.loc[
                        row['experiment_id'], 'data_path_cellmetadata'
                    ],
                    sep='\t',
                    index_col='cell_barcode'
                )

                for col in cellmetadata.columns:
                    adata.obs[col] = cellmetadata.loc[adata.obs.index, col]

        # Calculate basic qc metrics for this sample.
        # NOTE: n_genes_by_counts == number of genes with > 0 counts
        #       adata.obs['n_genes'] = (adata.X > 0).sum(axis = 1) is same as
        #       adata.obs['n_genes_by_counts']
        vars_prior_metrics = adata.var_keys()
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mito_gene'], inplace=True)

        # Apply cell filter.
        # adata = adata[selected_cells, :]
        # Apply gene filter
        # adata = adata[:, selected_genes]

        # Apply cell QC filters.
        adata.obs['cell_passes_qc'] = True
        filters_all_samples = []
        filters_experiment = []
        if 'cell_filters' not in params_dict:
            warnings.warn('Found no cell_filters in params_dict.')
        else:
            if 'all_samples' in params_dict['cell_filters'].keys():
                # NOTE: we want this to throw an error if value is not there.
                filters_all_samples = params_dict['cell_filters'][
                    'all_samples'
                ]['value']
            if row['experiment_id'] in params_dict['cell_filters'].keys():
                filters_experiment = params_dict['cell_filters'][
                    row['experiment_id']
                ]['value']

        # First record the total number of cells that pass each filter
        # independently i.e., not depenedent on any other filter.
        if len(filters_all_samples) > 0:
            for filter_query in filters_all_samples:
                if filter_query != '':
                    n_cells_dict[row['experiment_id']][
                        'filter__all_samples {}'.format(filter_query)
                    ] = adata.n_obs - adata.obs.query(filter_query).shape[0]
        if len(filters_experiment) > 0:
            for filter_query in filters_experiment:
                # NOTE: could add if test here to grab filter_query ==
                # file_cellids_filter or file_cellids_keep
                if filter_query != '':
                    n_cells_dict[row['experiment_id']][
                        'filter__sample_specific {}'.format(
                            filter_query
                        )
                    ] = adata.n_obs - adata.obs.query(filter_query).shape[0]

        # Now apply the filters - first apply the filters for all samples.
        n_cells_start = adata.n_obs
        filter_i = 0
        if len(filters_all_samples) > 0:
            # Run each filter iteratively.
            for filter_query in filters_all_samples:
                if filter_query != '':
                    # Drop the cells that are flagged in this query
                    cells_to_remove = adata.obs.query(filter_query).index
                    adata.obs.loc[cells_to_remove, 'cell_passes_qc'] = False
                    # adata = adata[
                    #     np.invert(adata.obs.index.isin(cells_to_remove)),
                    #     :
                    # ]
                    print('[{}] {} "{}": {} dropped {} remain'.format(
                        'all sample cell QC applied',
                        row['experiment_id'],
                        filter_query,
                        len(cells_to_remove),
                        adata.obs['cell_passes_qc'].sum()
                    ))
                    n_cells_dict[row['experiment_id']][
                        'filter__all_samples after_filter_{} {}'.format(
                            filter_i,
                            filter_query
                        )
                    ] = adata.obs['cell_passes_qc'].sum()
                    filter_i += 1

        # Now apply per sample filters.
        if len(filters_experiment) > 0:
            # Run each filter iteratively.
            for filter_query in filters_experiment:
                if filter_query != '':
                    cells_to_remove = adata.obs.query(filter_query).index
                    adata.obs.loc[cells_to_remove, 'cell_passes_qc'] = False
                    # adata = adata[
                    #     np.invert(adata.obs.index.isin(cells_to_remove)),
                    #     :
                    # ]
                    print('[{}] {} "{}": {} dropped {} remain'.format(
                        'sample specific cell QC applied',
                        row['experiment_id'],
                        filter_query,
                        len(cells_to_remove),
                        adata.obs['cell_passes_qc'].sum()
                    ))
                    n_cells_dict[row['experiment_id']][
                        'filter__sample_specific after_filter_{} {}'.format(
                            filter_i,
                            filter_query
                        )
                    ] = adata.obs['cell_passes_qc'].sum()
                    filter_i += 1

        # Write the number of cells filtered to standard out.
        print('[{}] after all cell QC: {} dropped {} remain'.format(
            row['experiment_id'],
            n_cells_start - adata.obs['cell_passes_qc'].sum(),
            adata.obs['cell_passes_qc'].sum()
        ))

        # Apply cell downsampling if needed.
        if params_dict['downsample_cells_fraction']['value'] != '':
            n_cells_start = adata.n_obs
            sc.pp.subsample(
                adata,
                fraction=float(
                    params_dict['downsample_cells_fraction']['value']
                ),
                copy=False,
                random_state=0
            )
            n_cells_dict[
                row['experiment_id']
            ]['downsample_cells_fraction'] = adata.n_obs
            print('[{}] cell downsample applied: {} dropped {} remain'.format(
                row['experiment_id'],
                n_cells_start - adata.n_obs,
                adata.n_obs
            ))
        elif params_dict['downsample_cells_n']['value'] != '':
            n_cells_start = adata.n_obs
            sc.pp.subsample(
                adata,
                n_obs=int(params_dict['downsample_cells_n']['value']),
                copy=False,
                random_state=0
            )
            n_cells_dict[
                row['experiment_id']
            ]['downsample_cells_n'] = adata.n_obs
            print('[{}] cell downsample applied: {} dropped {} remain'.format(
                row['experiment_id'],
                n_cells_start - adata.n_obs,
                adata.n_obs
            ))
        # Apply count downsampling if needed.
        if params_dict['downsample_feature_counts']['value'] != '':
            fraction = params_dict['downsample_feature_counts']['value']
            target_counts_per_cell = adata.obs['total_counts'].apply(
                lambda x: int(x * fraction)
            ).values
            sc.pp.downsample_counts(
                adata,
                counts_per_cell=target_counts_per_cell,
                random_state=0
            )

        # Print the number of cells and genes for this sample.
        n_cells_dict[row['experiment_id']]['after_filters'] = adata.obs[
            'cell_passes_qc'
        ].sum()
        print('[{}] {} obs (cells), {} var (genes)'.format(
            row['experiment_id'],
            adata.obs['cell_passes_qc'].sum(),
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

        # Only keep cells that pass QC
        if drop_filtered_cells:
            adata = adata[adata.obs['cell_passes_qc'], :]
            del adata.obs['cell_passes_qc']

        # If we still have cells after filters, add to our list of data.
        if adata.n_obs > 0:
            adatasets.append(adata)
            n_adatasets += 1
        else:
            raise Exception(
                'Error invalid tenx_data file. Missing coluns.'
            )

    # Merge all of the data together.
    adata_merged = adatasets[0].concatenate(*adatasets[1:])
    adata_merged = check_adata(adata_merged, 'adata_merged')

    # Re-calculate basic qc metrics of var (genes) for the whole dataset.
    # NOTE: we are only changing adata.var
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

    # Merge info on cell filters
    n_cells_df = pd.DataFrame(n_cells_dict)
    n_cells_df = n_cells_df.transpose()
    n_cells_df['experiment_id'] = n_cells_df.index
    n_cells_df = n_cells_df.melt(
        id_vars=['experiment_id'],
        var_name='filter_type',
        value_name='n_cells_left_in_adata'
    )
    # Drop rows with no value for n_cells_left_in_adata. This will happen for
    # per sample filters.
    n_cells_df = n_cells_df.dropna(subset=['n_cells_left_in_adata'])
    adata_merged.uns['cell_filtered_per_experiment'] = n_cells_df
    # adata_merged.uns['cell_filtered_per_experiment_dict'] = n_cells_dict

    # Save the adata matrix
    # output_file = output_dir + "/adata"
    adata_merged.write('{}.h5ad'.format(output_file), compression='gzip')
    # adata_merged.write_csvs(output_file)
    # adata_merged.write_loom(output_file+".loom")

    n_cells_df.to_csv(
        '{}-cell_filtered_per_experiment.tsv.gz'.format(output_file),
        sep='\t',
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        # index_label='cell_barcode',
        na_rep='',
        compression='gzip'
    )

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
        '-mf', '--sample_metadata_file',
        action='store',
        dest='mf',
        required=True,
        help='File with metadata on samples matching experiment_id column in\
            tenxdata_file. The column that links these two files is specified\
            in metadata_key.'
    )

    parser.add_argument(
        '-mcd', '--sample_metadata_columns_delete',
        action='store',
        dest='mcd',
        default='sample_status,study,study_id',
        help='Comma seperated list of columns to delete in\
            sample_metadata_file. If "" then no columns are deleted.\
            Note: whitespace should be represented with an underscore (_).\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-mk', '--metadata_key',
        action='store',
        dest='mk',
        default='experiment_id',
        help='Key to link metadata to tenxdata_file experiment_id column.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-cmf', '--cell_metadata_file',
        action='store',
        dest='cmf',
        default='',
        help='File with the following headers: experiment_id\
            data_path_cellmetadata, where data_path_cellmetadata is the path\
            to a tsv file contaning a cell_barcode column followed by other\
            columns to add to the annotations of each cell of that experiment\
            (e.g., doublet scores)\
            (default: None)'
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
            working directory. Will have .h5ad appended.\
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

    # Read in the parameters for downsampling and cell filters.
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

    # Load cell metadata
    cellmetadata_filepaths = None
    if options.cmf != '':
        cellmetadata_filepaths = pd.read_csv(
            options.cmf,
            sep='\t',
            index_col='experiment_id'
        )

    # Run the merge function.
    out_file = scanpy_merge(
        tenx_data,
        metadata,
        metadata_key=options.mk,
        output_file=options.of,
        params_dict=params_dict,
        cellmetadata_filepaths=cellmetadata_filepaths
    )
    print(out_file)


if __name__ == '__main__':
    main()
