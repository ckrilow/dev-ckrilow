#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
# import yaml

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
    samplesheetdata,
    metadata,
    output_file,
    sample_id_metadata_column="sanger_sample_id"
):
    """Merge 10x data.

    Parameters
    ----------
    samplesheetdata : pandas.DataFrame
        Description of parameter `samplesheetdata`.
    metadata : pandas.DataFrame
        Description of parameter `metadata`.
    output_file : string
        Description of parameter `output_file`.
    sample_id_metadata_column : string
        Description of parameter `sample_id_metadata_column`
        (the default is "sanger_sample_id").

    Returns
    -------
    output_file : string
        output_file
    """
    # check the samplesheetdata
    # check for required columns
    samplesheetdata_required_cols = set(['sample_id', 'data_path_10x_format'])
    if not samplesheetdata_required_cols.issubset(samplesheetdata.columns):
        raise Exception('Invalid samplesheetdata.')
    # check no duplicate sample ids
    vals, counts = np.unique(samplesheetdata['sample_id'], return_counts=True)
    if np.sum(counts > 1):
        raise Exception('Error {} duplicate sample_ids:\t{}'.format(
            np.sum(counts > 1),
            np.array2string(vals[counts > 1])
        ))
    # check all files exist
    filt = samplesheetdata['data_path_10x_format'].apply(
        lambda x: os.path.exists(x)
    )
    if np.sum(filt > 1):
        raise Exception('Error {} data_path_10x_format missing:\t{}'.format(
            np.sum(filt > 1),
            np.array2string(samplesheetdata['data_path_10x_format'][filt])
        ))

    # iterate over samples and load data
    adatasets = []
    n_adatasets = 1
    for idx, row in samplesheetdata.iterrows():
        # load the data
        adata = sc.read_10x_mtx(
            path=row['data_path_10x_format'],
            # var_names='gene_symbols',
            var_names='gene_ids',
            make_unique=False
        )

        adata = check_adata(adata, row['sample_id'])

        # label mitochondrial genes
        # mito_gene_list = sc.queries.mitochondrial_genes() # another query
        adata.var['mito_gene'] = [
            x.startswith('MT-') for x in adata.var['gene_symbols']
        ]
        # use this if var_names='gene_symbols' in sc.read_10x_mtx
        # adata.var['mito_gene'] = [
        #     x.startswith('MT-') for x in adata.var_names
        # ]

        # calculate basic qc metrics
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mito_gene'], inplace=True)

        # add in sample metadata
        # NOTE: it would be more memory efficient to stash this in
        #       unstructured dict-like annotation (adata.uns)
        metadata_smpl = metadata[
            metadata[sample_id_metadata_column] == row['sample_id']
        ]
        for col in metadata_smpl.columns:
            adata.obs[col] = np.repeat(metadata_smpl[col].values, adata.n_obs)

        # apply cell filter
        # adata = adata[selected_cells, :]
        # apply gene filter
        # adata = adata[:, selected_genes]

        # TODO
        # apply subsampling if specified in the config
        config_dict = dict()
        # config_dict['cell_qc'] = [
        #     'pct_counts_mito_gene < 80'
        # ]
        if 'subsample_cells' in config_dict.keys():
            sc.pp.subsample(adata, fraction=config_dict['subsample_cells'])
        if 'subsample_feature_counts' in config_dict.keys():
            fraction = config_dict['subsample_feature_counts']
            target_counts_per_cell = adata.obs['total_counts'].apply(
                lambda x: int(x * fraction)
            ).values
            sc.pp.downsample_counts(
                adata,
                counts_per_cell=target_counts_per_cell
            )
        if 'cell_qc' in config_dict.keys():
            n_cells_start = adata.n_obs
            for filter_query in config_dict['cell_qc']:
                adata = adata[adata.obs.query(filter_query).index, :]
                print('[{}] cell QC applied "{}": {} cells dropped'.format(
                    row['sample_id'],
                    filter_query,
                    n_cells_start - adata.n_obs
                ))
            print('[{}] cell QC applied: {} total cells dropped'.format(
                row['sample_id'],
                n_cells_start - adata.n_obs
            ))

        # print the number of cells and genes for this sample
        print('[{}] {} obs (cells), {} vars (genes)'.format(
            row['sample_id'],
            adata.n_obs,
            adata.n_vars
        ))

        # if we still have cells after filters, add to our list of data
        if adata.n_obs > 0:
            adatasets.append(adata)
            n_adatasets += 1

    # merge all of the data together
    adata_merged = adatasets[0].concatenate(*adatasets[1:])
    adata_merged = check_adata(adata_merged, 'adata_merged')

    # possible additional basic filtering on the full dataset
    # sc.pp.filter_cells(adata, min_genes=200)
    # sc.pp.filter_genes(adata, min_cells=3)

    print('[adata_merged] {} obs, {} vars'.format(
        adata_merged.n_obs,
        adata_merged.n_vars
    ))

    # output_file = output_dir + "/adata"
    adata_merged.write(output_file+'.h5', compression='gzip')
    # adata_merged.write_csvs(output_file)
    # adata_merged.write_loom(output_file+".loom")

    return(output_file)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter and merge 10x data. Save to anndata object.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-sf', '--samplesheetdata_file',
        action='store',
        dest='sf',
        required=True,
        help='File with the following headers: sanger_sample_id, file_path.'
        # 'Input yaml config file that lists the samples to merge as well',
        # ' as filters to apply to each file.'
    )

    parser.add_argument(
        '-mf', '--metadata_file',
        action='store',
        dest='mf',
        required=True,
        help='File with metadata on samples matching sanger_sample_id in\
            samplesheetdata_file.'
    )

    parser.add_argument(
        '-od', '--output_dir',
        action='store',
        dest='od',
        default=os.getcwd(),
        help='Directory to write output anndata file.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # requirements:
    # - reads yaml config file with a list per sample and the filters one
    #   wants to use for that sample
    # - input file can be either 10x dir or a loom object. Use loom object
    #   to add metadata to samples prior to merge that might be used in
    #   the filtering.
    # output:
    # - each sample file is read in, filtered, and added to a massive dataframe
    #   at the end the dataframe is written out - Q: possible to not hold all
    #   dataset in memory by writing to loom file one at a time?
    # with open(config_file, 'r') as f:
    #     config_dict = yaml.safe_load(f)

    # load a file of the samples to analyse
    # sf = "/home/ubuntu/studies/TaylorDL/freeze001/final_samples.tsv"
    samplesheetdata = pd.read_csv(options.sf, sep='\t')
    # samplesheetdata = samplesheetdata.drop(['experiment_id'], axis=1)
    samplesheetdata = samplesheetdata.rename(columns={
        'sanger_sample_id': 'sample_id',
        'file_path': 'data_path_10x_format'
    })

    # load the metadata
    # f = "/home/ubuntu/repo/scrna_cellranger/sync_status/samples_metainfo.tsv"
    metadata_column_delete = ['sample_status', 'study', 'study_id']
    metadata = pd.read_csv(options.mf, sep='\t')
    metadata.columns = metadata.columns.str.strip(
        ).str.replace(' ', '_').str.lower()
    metadata = metadata.drop(metadata_column_delete, axis=1)

    out_file = scanpy_merge(
        samplesheetdata,
        metadata,
        output_file=options.od.rstrip('/') + '/adata'
    )
    print(out_file)


if __name__ == '__main__':
    main()
