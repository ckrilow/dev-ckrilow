#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import csv
import time
from datetime import timedelta

# Set scanpy settings
# sc verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.settings.verbosity = 3
# sc.logging.print_versions()
# sc.settings.set_figure_params(dpi=80)


def score_cells(
    adata,
    score_genes_df,
    score_genes_df_column='ensembl_gene_id',
    only_use_variable_genes=False
):
    """Scores each cell.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object. Assume adata.X is norm->log1p->scaled data.
    score_genes_df : pd.DataFrame
        Dataframe of marker genes. Needs to have score_genes_df_column and
        score_id column. If one score_id == 'cell_cycle', then requires a
        grouping_id column with 'G2/M' and 'S'.
    score_genes_df_column : string
        Column in score_genes_df to use for gene ids (e.g., hgnc_symbol,
        ensembl_gene_id)
    only_use_variable_genes : boolean
        Only use variable genes to calculate scores. If True, score_id will
        be changed to <score_id>__hvg_only. Note this flage does not apply
        to score_id == 'cell_cycle'.


    Returns
    -------
    adata : AnnData
        AnnData object with scores calculated and stored in
        adata.obs[<score_id>].
    score_genes_df : pd.DataFrame
        The score_genes_df with the following columns added:
        gene_found_in_adata, gene_found_is_highly_variable. It is suggested
        that this dataframe is added to the adata.uns slot.
    """
    verbose = False  # For debugging purposes.

    # Update the score_genes_df with details on the genes and if they were
    # found in adata and if they are highly variable.
    score_genes_df['gene_found_in_adata'] = np.in1d(
        score_genes_df[score_genes_df_column],
        adata.var.index
    )
    score_genes_df['gene_found_is_highly_variable'] = np.in1d(
        score_genes_df[score_genes_df_column],
        adata.var.index[adata.var['highly_variable']]
    )

    # Set the gene pool parameter.
    gene_pool = None  # If None, all genes are randomly sampled for background
    if only_use_variable_genes:
        gene_pool = adata.var.index[adata.var['highly_variable']]

    # Loop over each score_id in score_genes_df, updating adata.
    for score_id, df_group in score_genes_df.groupby('score_id'):
        # Downsample to only those genes found in the data.
        df_group = df_group.loc[
            df_group['gene_found_in_adata'], :
        ]
        if df_group.shape[0] == 0:
            continue

        # If we are supposed to use only_use_variable_genes, then do so.
        if only_use_variable_genes:
            if score_id == 'cell_cycle':
                continue
            score_id = '{}__hvg_only'.format(score_id)
            df_group = df_group.loc[
                df_group['gene_found_is_highly_variable'], :
            ]
            if df_group.shape[0] == 0:
                continue
        if verbose:
            print('Scoring {}'.format(score_id))

        # Set the number of control genes.
        ctrl_size = 50
        if df_group.shape[0] > 50:
            ctrl_size = df_group.shape[0]
        if gene_pool is not None:
            if ctrl_size > len(gene_pool):
                raise Exception(
                    'Error in gene scoring ctrl_size > len(gene_pool)'
                )

        # If the score_id is cell_cycle, then use the specific cell cycle
        # scoring function.
        if score_id == 'cell_cycle':
            # NOTE: Setting ctrl_size` is not possible, as it's set as
            #       `min(len(s_genes), len(g2m_genes))`.
            sc.tl.score_genes_cell_cycle(
                adata,
                s_genes=df_group.loc[
                    df_group['grouping_id'] == 'S', score_genes_df_column
                ],
                g2m_genes=df_group.loc[
                    df_group['grouping_id'] == 'G2/M', score_genes_df_column
                ],
                copy=False,
                gene_pool=gene_pool,  # Default is None (aka, use all)
                n_bins=25,  # Default is 25
                use_raw=False
            )
        else:
            sc.tl.score_genes(
                adata,
                df_group[score_genes_df_column],
                ctrl_size=ctrl_size,  # Default is 50
                gene_pool=gene_pool,  # Default is None (aka, use all)
                n_bins=25,  # Default is 25
                score_name=score_id,
                random_state=0,  # Default is 0
                copy=False,
                use_raw=False
            )

    return adata, score_genes_df


def scanpy_normalize_and_pca(
    adata,
    output_file,
    vars_to_regress,
    variable_feature_batch_key='sanger_sample_id',
    n_variable_features=2000,
    exclude_hv_gene_df=[],
    score_genes_df=None,
    verbose=True,
    plot=True
):
    """Normalize data and calculate PCs.

    Parameters
    ----------
    adata : AnnData
        Input AnnData file.
    output_file : string
        Basename of output_file, will have -normalized_pca.h5ad appended to it.
    vars_to_regress : list
        List of metadata variables to regress. If empty no regression.
    variable_feature_batch_key : string
        Batch key for variable gene detection.
        The default is "sanger_sample_id".
    n_variable_features : int
        Number of variable features to select.
    exclude_hv_gene_df : pd.DataFrame
        Dataframe of genes to exclude from highly variable gene selection.
    score_genes_df : pd.DataFrame
        Dataframe of marker genes. Needs to have score_genes_df_column and
        score_id column. If one score_id == 'cell_cycle', then requires a
        grouping_id column with 'G2/M' and 'S'.
    verbose : boolean
        Write extra info to standard out.
    plot : boolean
        Generate plots.


    Returns
    -------
    output_file : string
        output_file
    """
    # Check that any vars to regress occur in adata
    if len(vars_to_regress) > 0:
        for i in vars_to_regress:
            if i not in adata.obs.columns:
                raise Exception(
                    '{} in vars_to_regress missing from metadata'.format(
                        i
                    )
                )

    # Add a raw counts layer.
    # NOTE: This stays with the main AnnData and is not stashed when we
    #       later save the ln(CPM+1) data to raw (raw only stores X without
    #       layers).
    adata.layers['counts'] = adata.X.copy()

    # NOTE: prior to running normalization, low quality cells should be
    # filtered. Example:
    # sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=5)
    # Only consider genes expressed in more than 0.5% of cells:
    # sc.pp.filter_genes(adata, min_cells=0.005*len(adata.obs.index))

    # Total-count normalize (library-size correct) the data matrix X to
    # counts per million, so that counts become comparable among cells.
    sc.pp.normalize_total(
        adata,
        target_sum=1e6,
        exclude_highly_expressed=False,
        key_added='normalization_factor',  # add to adata.obs
        inplace=True
    )
    # Logarithmize the data: X = log(X + 1) where log = natural logorithm.
    # Numpy has a nice function to undo this np.expm1(adata.X).
    sc.pp.log1p(adata)
    # Delete automatically added uns
    del adata.uns['log1p']
    # Add record of this operation.
    adata.layers['log1p_cpm'] = adata.X.copy()
    adata.uns['log1p_cpm'] = {'transformation': 'ln(CPM+1)'}

    # Stash the unprocessed data in the raw slot.
    # adata.raw.X.data is now ln(CPM+1).
    # NOTE: - Layers are not preserved in adata.raw, though obs, var, uns are.
    #       - If genes are filtered (e.g.,
    #         sc.pp.filter_genes(adata, min_cells=1)), the full dataset will
    #         remain in the raw slot.
    #       - We store in the raw slot because later for UMAP and marker gene
    #         analysis, we can easily tell scanpy to use the raw slot via the
    #         use_raw = True flag. Raw was specifically designed for this use
    #         case of ln(CPM+1),
    # Can be deleted later: del adata.raw
    adata.raw = adata
    # adata_raw = adata.raw.to_adata()

    if plot:
        # Plot top expressed genes.
        _ = sc.pl.highest_expr_genes(
            # adata.raw.to_adata(),  # same as adata at this point.
            adata,
            n_top=25,
            gene_symbols='gene_symbols',
            show=False,
            save='-{}.pdf'.format(output_file)
        )

    # Calculate the highly variable genes on the log1p(norm) data.
    # Do so for each sample and then merge - this avoids the selection of
    # batch-specific, highly variable genes.
    sc.pp.highly_variable_genes(
        adata,
        # min_mean=0.0125,
        # max_mean=3,
        # min_disp=0.5,
        flavor='seurat',
        n_top_genes=n_variable_features,  # 2000 = SeuratFindVariableFeatures
        batch_key=variable_feature_batch_key,
        inplace=True
    )
    if verbose:
        print('{}: {} (all batches); {} ({})'.format(
            'Number of variable features detected',
            adata.var['highly_variable_intersection'].sum(),
            adata.var['highly_variable'].sum(),
            'after ranking the number of batches where a feature is variable'
        ))
    # If n_top_genes = None, then one needs to set 'highly_variable'.
    # Here, highly_variable_intersection is only true for genes variable across
    # all batch keys (i.e., 'highly_variable_nbatches' = n_batch_keys):
    # adata.var.loc[
    #     adata.var["highly_variable_intersection"],
    #     ["highly_variable_nbatches"]
    # ]
    #
    # If n_top_genes = None, then one also needs needs to set highly_variable'.
    # Fix bug in PCA when we have set batch_key. More below:
    # https://github.com/theislab/scanpy/issues/1032
    # adata.var['highly_variable'] = adata.var['highly_variable_intersection']
    #
    # Alternatively, if one specifies n_top_genes, then genes are ranked by
    # 'highly_variable_nbatches' and highly_variable is set to the top n.
    # adata.var.loc[
    #     adata.var["highly_variable"],
    #     ["highly_variable_nbatches"]
    # ]

    if plot:
        # Plot highly variable genes.
        _ = sc.pl.highly_variable_genes(
            adata,
            log=False,
            show=False,
            save='-{}.pdf'.format(output_file)
        )
        _ = sc.pl.highly_variable_genes(
            adata,
            log=True,
            show=False,
            save='-{}-log.pdf'.format(output_file)
        )

    # After calculating highly variable genes, we subsquently remove any custom
    # for highly variable gene selection. This way we retain the normalized
    # values for each one of these genes even though they will not be used
    # for dimensionality reduction. NOTE: If there are loads of genes
    # to exclude and there are only a handful of n_variable_features, then
    # one could end up with very few variable genes for dimensionality
    # reduction in the end.
    #
    # Exclude mitocondrial genes from highly variable gene set.
    # if exclude_mito_highly_variable_genes:
    #     n_highly_variable_mito = adata.var.loc[
    #         adata.var['mito_gene'],
    #         ['highly_variable']
    #     ].sum()
    #     if verbose:
    #         print('Within highly variable genes, {} are mito genes'.format(
    #             n_highly_variable_mito
    #         ))
    #     adata.var.loc[
    #         adata.var['mito_gene'],
    #         ['highly_variable']
    #     ] = False
    # Exclude other genes from highly variable gene set.
    if len(exclude_hv_gene_df) > 0:
        # Annotate the exclusion dataframe with the genes that are highly
        # variable.
        exclude_hv_gene_df['highly_variable'] = exclude_hv_gene_df[
            'ensembl_gene_id'
        ].isin(adata.var.loc[adata.var.highly_variable, :].index)

        # Exclude these genes.
        adata.var.loc[
            exclude_hv_gene_df.loc[
                exclude_hv_gene_df.highly_variable, :
            ]['ensembl_gene_id'],
            ['highly_variable']
        ] = False

        # Add record of gene exclustions
        adata.uns['df_highly_variable_gene_filter'] = exclude_hv_gene_df

        # Print out the number of genes excluded
        if verbose:
            print('Within highly variable genes, {} genes are {}'.format(
                exclude_hv_gene_df['highly_variable'].sum(),
                'in the list of genes to exclude.'
            ))

    if len(vars_to_regress) == 0:
        # Scale the data to unit variance.
        # This effectively weights each gene evenly.
        sc.pp.scale(
            adata,
            zero_center=True,
            max_value=None,
            copy=False
        )
        # Calculate gene scores on each cell.
        # Perform this two ways:
        # (1) All genes [that pass basic 0 filters]. As done in sc tutorial:
        # https://github.com/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
        # (2) Only highly variable genes. As done in:
        # https://www.biorxiv.org/content/10.1101/2020.04.03.024075v1
        if score_genes_df is not None:
            adata, score_genes_df_updated = score_cells(
                adata,
                score_genes_df,
                score_genes_df_column='ensembl_gene_id',
                only_use_variable_genes=False
            )
            adata, _ = score_cells(
                adata,
                score_genes_df,
                score_genes_df_column='ensembl_gene_id',
                only_use_variable_genes=True
            )
    else:  # Regress out any continuous variables.
        # Before regressing calculate the gene scores on a copy of the data.
        if score_genes_df is not None:
            adata_scored = sc.pp.scale(
                adata,
                zero_center=True,
                max_value=None,
                copy=True
            )
            adata_scored, score_genes_df_updated = score_cells(
                adata_scored,
                score_genes_df,
                score_genes_df_column='ensembl_gene_id',
                only_use_variable_genes=False
            )
            adata_scored, _ = score_cells(
                adata_scored,
                score_genes_df,
                score_genes_df_column='ensembl_gene_id',
                only_use_variable_genes=True
            )

            # Add scores back into the main dataframe.
            new_cols = np.setdiff1d(
                adata_scored.obs.columns,
                adata.obs.columns
            )
            adata.obs = pd.concat(
                [adata.obs, adata_scored.obs.loc[adata.obs.index, new_cols]],
                axis=1
            )

        # NOTE: if the same value is repeated (e.g., 0) for all cells this will
        #       fail. https://github.com/theislab/scanpy/issues/230
        # if verbose:
        #     print('For regress_out, calling {}'.format(
        #         'pp.filter_genes(adata, min_cells=5)'
        #     ))
        # sc.pp.filter_genes(adata, min_cells=5)
        # NOTE: sc.pp.regress_out out should default to sc.settings.n_jobs
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
            max_value=None,
            copy=False
        )

    # Keep a record of the different gene scores
    if score_genes_df is not None:
        adata.uns['df_score_genes'] = score_genes_df_updated

    # Calculate PCs.
    sc.tl.pca(
        adata,
        n_comps=min(200, adata.var['highly_variable'].sum()),
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
        '{}-pcs.tsv.gz'.format(output_file),
        sep='\t',
        index=True,
        index_label='cell_barcode',
        na_rep='',
        compression='gzip'
    )

    # Save the metadata to a seperate file for Harmony.
    adata.obs.to_csv(
        '{}-metadata.tsv.gz'.format(output_file),
        sep='\t',
        index=True,
        quoting=csv.QUOTE_NONNUMERIC,
        index_label='cell_barcode',
        na_rep='',
        compression='gzip'
    )

    # Save the data.
    adata.write(
        '{}-normalized_pca.h5ad'.format(output_file),
        compression='gzip'
    )
    # adata_merged.write_csvs(output_file)
    # adata_merged.write_loom(output_file+".loom"))

    # Plot the PC info.
    if plot:
        # Plot the vanilla PCs.
        # sc.pl.pca(
        #     adata,
        #     color='sanger_sample_id',
        #     components=['1,2', '3,4']
        # )
        _ = sc.pl.pca_variance_ratio(
            adata,
            n_pcs=adata.obsm['X_pca'].shape[1],
            log=False,
            show=False,
            save='-{}.pdf'.format(output_file)
        )
        _ = sc.pl.pca_variance_ratio(
            adata,
            n_pcs=adata.obsm['X_pca'].shape[1],
            log=True,
            show=False,
            save='-{}-log.pdf'.format(output_file)
        )

    # Save the filtered count matrix for input to other software like scVI
    adata.X = adata.layers['counts']
    del adata.layers['counts']
    del adata.raw
    adata.write(
        '{}-normalized_pca-counts.h5ad'.format(output_file),
        compression='gzip'
    )

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
        '-bk', '--batch_key',
        action='store',
        dest='bk',
        default='sanger_sample_id',
        help='Batch key for highly-variable feature (e.g., gene) detection.\
            If specified, highly-variable features are selected within each\
            batch separately and merged.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-nvf', '--number_variable_features',
        action='store',
        dest='nvf',
        default=2000,
        type=int,
        help='After calculating variable features within each batch set via\
            <batch_key>, rank features by number of batches where they are\
            variable and select the top <number_variable_features>.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-vge', '--variable_genes_exclude',
        action='store',
        dest='vge',
        default='',
        help='Tab-delimited file with genes to exclude from the highly\
            variable gene list. Must contain ensembl_gene_id column.\
            (default: None - keep all variable genes)'
    )

    parser.add_argument(
        '-vr', '--vars_to_regress',
        action='store',
        dest='vr',
        default='',
        help='Comma seperated list of metadata variables to regress prior to\
            calculating PCs. Example: mito_gene,n_count.\
            (default: "" and sc.pp.regress_out is not called)'
    )

    parser.add_argument(
        '-sg', '--score_genes',
        action='store',
        dest='sg',
        default='',
        help='Tab-delimited file of genes for scores. Needs to have\
            ensembl_gene_id and score_id column. If one\
            score_id == "cell_cycle", then requires a grouping_id column with\
            "G2/M" and "S".'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=4,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='adata-normalize_pca',
        help='Directory and basename of output files.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    # Split the vars to regress list
    vars_to_regress = []
    if options.vr != '':
        vars_to_regress = options.vr.split(',')

    # Load list of genes to filter
    genes_filter = []
    if options.vge != '':
        genes_filter = pd.read_csv(options.vge, sep='\t')

    # Load the gene scores
    score_genes_df = None
    if options.sg != '':
        score_genes_df = pd.read_csv(options.sg, sep='\t')

    start_time = time.time()
    _ = scanpy_normalize_and_pca(
        adata,
        output_file=options.of,
        vars_to_regress=vars_to_regress,
        variable_feature_batch_key=options.bk,
        n_variable_features=options.nvf,
        exclude_hv_gene_df=genes_filter,
        score_genes_df=score_genes_df,
        verbose=True
    )
    execution_summary = "Analysis execution time [{}]:\t{}".format(
        "scanpy_normalize_and_pca.py",
        str(timedelta(seconds=time.time()-start_time))
    )
    print(execution_summary)


if __name__ == '__main__':
    main()
