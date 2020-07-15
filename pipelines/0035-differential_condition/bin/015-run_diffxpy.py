#!/usr/bin/env python


__author__ = 'Monika Krzak and Leland Taylor'
__date__ = '2020-07-02'
__version__ = '0.0.1'

import argparse
import warnings
import numpy as np
import scanpy as sc
import diffxpy.api as de
import joblib  # for numpy matrix, joblib faster than pickle


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Differential expression analysis using Wald test.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-cond', '--condition_column',
        action='store',
        dest='condition_column',
        required=True,
        help='Condition column. Should be a binary trait.'
    )

    parser.add_argument(
        '-cov', '--covariate_columns',
        action='store',
        dest='cov',
        default='normalization_factor',
        help='Comma seperated list of covariates (columns in adata.obs). \
            (default: %(default)s)'
    )

    parser.add_argument(
        '-cl', '--cell_label_column',
        action='store',
        dest='cell_label_column',
        default='cluster',
        help='Anndata cell type label name in obs slot. (default: %(default)s)'
    )

    parser.add_argument(
        '-cla', '--cell_label_analyse',
        action='store',
        dest='cell_label_analyse',
        default='',
        help='Only analyse cells with these labels. Either "all" or comma \
            seperated list of cell labels. (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output file.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get the parameters
    condition_column = options.condition_column
    covariate_columns = options.covariate_columns
    cell_label_column = options.cell_label_column
    cell_label_analyse = options.cell_label_analyse
    out_file_base = options.of

    # Diffxpy automatically assigns covariates to categorical. List to stash
    # continuous variables.
    continuous_variables = []

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)
    adata.X = adata.layers['counts']

    # Sort out the condition we want to test
    adata.obs['condition'] = adata.obs[condition_column]
    # Check to see if the condition is a categorical variable
    if adata.obs[condition_column].dtype.name != 'category':
        # raise Exception('Condition is not a category.')
        continuous_variables.append('condition')
        warnings.warn('Treating condition as a continuous variable')
    # The below code converts the condition into an int and asserts it binary
    # adata.obs['condition'] = adata.obs[condition_column].cat.codes
    # if len(np.unique(adata.obs['condition']) != 2):
    #     raise Exception('There are not exactly 2 conditions.')

    # Calculate size factor covariate manually.
    # NOTE: This is commented out and instead we use the 'normalization_factor'
    #       slot in the anndata.
    #
    # When running diffxpy regress size_factors. For more details check:
    # https://nbviewer.jupyter.org/github/theislab/diffxpy_tutorials/blob/master
    # /diffxpy_tutorials/test/introduction_differential_testing.ipynb
    # "Inclusion of continuous effects"
    # adata.obs['size_factors'] = adata.obs['total_counts']/10000

    # Figure out continuous covariates.
    covariate_columns = covariate_columns.split(',')
    for i in covariate_columns:
        if adata.obs[i].dtype.name != 'category':
            continuous_variables.append(i)

    # Select the subset of cells we want for analysis
    if cell_label_analyse != 'all':
        cell_label_analyse = cell_label_analyse.split(',')
        adata = adata[adata.obs[cell_label_column].isin(cell_label_analyse)]
    # clusters = np.sort(adata.obs[cell_label_column].unique())

    # Check to make sure that within each cluster,
    n_cells_condition_cluster = adata.obs.groupby(
        [condition_column, cell_label_column]
    ).size().unstack()
    if n_cells_condition_cluster.values.min() == 0:
        raise Exception('For one cell_label_column there are 0 conditions.')
    if len(np.unique(adata.obs['condition']) <= 1):
        raise Exception('There is only 1 condition.')

    print('Continuous varibles: {}'.format(','.join(continuous_variables)))

    # Run diffxpy
    #
    # Tutorial for conditions:
    # https://nbviewer.jupyter.org/github/theislab/diffxpy_tutorials/blob/master/diffxpy_tutorials/test/introduction_differential_testing.ipynb
    # Tutorial for continuous:
    # https://nbviewer.jupyter.org/github/theislab/diffxpy_tutorials/blob/master/diffxpy_tutorials/test/modelling_continuous_covariates.ipynb
    formula = '~ 1 + condition'
    for cov in covariate_columns:
        formula = '{} + {}'.format(formula, cov)
    # Get the coefficient names.
    # Example if coefficient = M or F then...
    # coef_to_test=['condition[T.M]']
    # If T1, T2, T3
    # coef_to_test=['condition[T.1]', 'condition[T.2]', 'condition[T.3]']
    coef_names = de.utils.preview_coef_names(
        sample_description=adata.obs,
        formula=formula,
        as_numeric=continuous_variables
    )
    coef_names = [i for i in coef_names if i.startswith('condition')]
    # param as_numeric:
    # Which columns of sample_description to treat as numeric and
    # not as categorical. This yields columns in the design matrix
    # which do not correspond to one-hot encoded discrete factors.
    # This makes sense for number of genes, time, pseudotime or space
    # for example.
    test_results = de.test.wald(
        data=adata,
        formula_loc=formula,
        factor_loc_totest=coef_names,
        as_numeric=continuous_variables
    )

    # Now save the results
    out_f = '{}-de_model.joblib.gz'.format(out_file_base)
    joblib.dump(
        test_results,
        out_f,
        compress=('gzip', 3)
    )

    # NOTE: the below code partitions the data within diffxpy ... we don't
    # need to do this since we have already partitioned clusters via the
    # cell_label_analyse command
    # part = de.test.partition(data=adata, parts="cluster")
    #
    # if options.cov == "":
    #     formula = "~ 1 + condition"
    #     test_part = part.wald(formula_loc=formula,
    #                           factor_loc_totest="condition")
    # else:
    #     formula = "~ 1 + condition + " + options.cov
    #     test_part = part.wald(formula_loc=formula,
    #                           factor_loc_totest="condition",
    #                           as_numeric=[options.cov])
    # for i in clusters:
    #     df = test_part.tests[test_part.partitions.index(i)].summary()
    #     df.to_csv('{}-cluster_de_wald.tsv'.format(i), sep='\t', index=False)
    #
    # part = de.test.partition(data=adata, parts="cluster")
    # test_part = part.wald(formula_loc="~ 1 + condition",
    #                       factor_loc_totest="condition")

    # BEGIN: developmental code
    # for i in clusters:
    #     df = test_part.tests[test_part.partitions.index(i)].summary()
    #     df.to_csv('{}-cluster_de_wald.tsv'.format(i), sep='\t', index=False)
    #
    # size_factors = adata.obs['total_counts']/10000
    # adata.obs['size_factors'] = size_factors
    # clusters = clusters[35:37]
    # for i in clusters:
    #
    #     adata1 = adata[adata.obs['cluster'] == i]
    #     part = de.test.partition(data=adata1, parts="cluster")
    #     test_part = part.wald(formula_loc="~ 1 + condition + size_factors",
    #                           factor_loc_totest="condition",
    #                           as_numeric=["size_factors"])
    #
    #     df = test_part.tests[test_part.partitions.index(i)].summary()
    #     df.to_csv('{}-cluster_de_wald.tsv'.format(i), sep='\t', index=False)
    #
    # size_factors = adata.obs['total_counts']/10000
    # adata.obs['size_factors'] = size_factors
    # part = de.test.partition(data=adata, parts="cluster")
    # test_part = part.wald(formula_loc="~ 1 + condition + size_factors",
    #                       factor_loc_totest="condition",
    #                       as_numeric=["size_factors"])
    # for i in clusters:
    #     df = test_part.tests[test_part.partitions.index(i)].summary()
    #     df.to_csv('{}-cluster_de_wald.tsv'.format(i), sep='\t', index=False)


if __name__ == '__main__':
    main()
