#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import pandas as pd
import scanpy as sc

# Set scanpy settings
# sc verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.settings.verbosity = 3
# sc.logging.print_versions()
# sc.settings.set_figure_params(dpi=80)


def scanpy_normalize_and_pca(
    adata,
    output_file,
    vars_to_regress,
    sample_id_metadata_column='sanger_sample_id',
    verbose=True
):
    """Normalize data and calculate PCs.

    Parameters
    ----------
    adata : AnnData
        Description of parameter `adata`.
    output_file : string
        Description of parameter `output_file`.
    vars_to_regress : list
        Description of parameter `vars_to_regress`.
    sample_id_metadata_column : string
        Description of parameter `sample_id_metadata_column`
        (the default is "sanger_sample_id").
    verbose : boolean
        Write extra info to standard out.

    Returns
    -------
    output_file : string
        output_file
    """
    # Stash the unprocessed data to the raw slot. Can be deleted later:
    # del adata.raw
    adata.raw = adata

    # NOTE: prior to running normalization, low quality cells should be
    # filtered. Example:
    # sc.pp.filter_cells(adata, min_genes=200)
    # sc.pp.filter_genes(adata, min_cells=3)
    # Only consider genes expressed in more than 0.5% of cells:
    # sc.pp.filter_genes(adata, min_cells=0.005*len(adata.obs.index))

    # Total-count normalize (library-size correct) the data matrix X to 10,000
    # reads per cell, so that counts become comparable among cells.
    sc.pp.normalize_total(
        adata,
        target_sum=1e4,
        exclude_highly_expressed=False,
        inplace=True
    )
    # Logarithmize the data: X = log(X + 1) where log = natural logorithm
    sc.pp.log1p(adata)

    # adata.raw.X.data is still raw count data
    # adata.X.data is now log normalized data

    # Calculate the highly variable genes. Do so for each sample and then
    # merge - this avoids the selection of batch-specific genes.
    sc.pp.highly_variable_genes(
        adata,
        # min_mean=0.0125,
        # max_mean=3,
        # min_disp=0.5,
        flavor='seurat',
        batch_key=sample_id_metadata_column,
        inplace=True
    )
    if verbose:
        print('number of variable genes detected across batches:\t{}'.format(
            adata.var['highly_variable_intersection'].sum()
        ))
    # Fix bug in PCA when we have set batch_key. More below:
    # https://github.com/theislab/scanpy/issues/1032
    adata.var['highly_variable'] = adata.var['highly_variable_intersection']

    # Tegress out any continuous variables.
    if (len(vars_to_regress) > 0):
        sc.pp.regress_out(
            adata,
            keys=vars_to_regress,
            copy=False
        )

    # Scale the data to unit variance.
    # This effectively weights each gene evenly.
    sc.pp.scale(
        adata,
        zero_center=True,
        max_value=None
    )

    # calculate PCs
    sc.tl.pca(
        adata,
        n_comps=200,
        zero_center=None,
        svd_solver='auto',
        use_highly_variable=True,
        copy=False
    )

    # Save PCs to a seperate file for Harmony.
    pca_df = pd.DataFrame(
        adata.obsm['X_pca'],
        index=adata.obs_names,
        columns=[
            'PC{}'.format(x) for x in range(1, adata.obsm['X_pca'].shape[1]+1)
        ]
    )
    pca_df.to_csv(
        output_file+'-pcs.csv.gz',
        sep=',',
        index=True,
        index_label='cell_barcode',
        compression='gzip'
    )

    # Save the metadata to a seperate file for Harmony.
    adata.obs.to_csv(
        output_file+'-metadata.csv.gz',
        sep=',',
        index=True,
        index_label='cell_barcode',
        compression='gzip'
    )

    # Plot top expressed genes.
    # sc.pl.highest_expr_genes(adata, n_top=20, gene_symbols='gene_symbols')
    # Plot highly variable genes.
    # sc.pl.highly_variable_genes(adata, log=False)
    # Plot the vanilla PCs.
    # sc.pl.pca(
    #     adata,
    #     color='sanger_sample_id',
    #     components=['1,2', '3,4']
    # )
    # sc.pl.pca_variance_ratio(adata, log=False)

    # Save the data.
    adata.write(output_file+'-normalized_pca.h5', compression='gzip')
    # adata_merged.write_csvs(output_file)
    # adata_merged.write_loom(output_file+".loom")

    return(output_file)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read anndata object. Normalize, calculate PCs. Save new anndata
            object along with csv file of PCs.
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
        '-vr', '--vars_to_regress',
        action='store',
        dest='vr',
        default='',
        help='Comma seperated list of metadata variables to regress prior to\
            calculating PCs. Example: mito_gene,n_count\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default=os.getcwd()+'/adata',
        help='Directory and basename of output files.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    # Split the vars to regress list
    vars_to_regress = []
    if options.vr != '':
        vars_to_regress = options.vr.split(',')

    out_file = scanpy_normalize_and_pca(
        adata,
        output_file=options.of,
        vars_to_regress=vars_to_regress,
        sample_id_metadata_column='sanger_sample_id',
        verbose=True
    )
    print(out_file)


if __name__ == '__main__':
    main()
