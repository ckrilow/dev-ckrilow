#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-07-02'
__version__ = '0.0.1'

import argparse
import numpy as np
import scanpy as sc
import diffxpy.api as de


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


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
            '-cov', '--covariate',
            action='store',
            dest='cov',
            required=True,
            help='Covariate to be regressed out. Can be "" or\
            name of the column that is in adata.obs'
    )

    parser.add_argument(
            '-ic', '--include_clusters',
            action='store',
            dest='ic',
            required=True,
            help='Clusters to be included. \
            Specify "all" or cluster number'
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

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)
    condition = np.transpose(
      np.array([adata.obs['disease_status'] == "Crohn's disease"]).astype(int)
    )
    adata.obs['condition'] = condition
    adata.X = adata.layers['counts'].copy()
    clusters = sorted(adata.obs['cluster'].unique(), key=int)

    # When running diffxpy regress size_factors. For more details check:
    # https://nbviewer.jupyter.org/github/theislab/diffxpy_tutorials/blob/master
    # /diffxpy_tutorials/test/introduction_differential_testing.ipynb
    # "Inclusion of continuous effects"

    size_factors = adata.obs['total_counts']/10000
    adata.obs['size_factors'] = size_factors

    # part = de.test.partition(data=adata, parts="cluster")
    # test_part = part.wald(formula_loc="~ 1 + condition",
    #                       factor_loc_totest="condition")
    #
    # for i in clusters:
    #     df = test_part.tests[test_part.partitions.index(i)].summary()
    #     df.to_csv('{}-cluster_de_wald.tsv'.format(i), sep='\t', index=False)

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

    # size_factors = adata.obs['total_counts']/10000
    # adata.obs['size_factors'] = size_factors
    # part = de.test.partition(data=adata, parts="cluster")
    # test_part = part.wald(formula_loc="~ 1 + condition + size_factors",
    #                       factor_loc_totest="condition",
    #                       as_numeric=["size_factors"])
    # for i in clusters:
    #     df = test_part.tests[test_part.partitions.index(i)].summary()
    #     df.to_csv('{}-cluster_de_wald.tsv'.format(i), sep='\t', index=False)

    if options.ic != "all":
        adata = adata[adata.obs['cluster'] == options.ic]
        clusters = sorted(adata.obs['cluster'].unique(), key=int)

    part = de.test.partition(data=adata, parts="cluster")

    if options.cov == "":
        formula = "~ 1 + condition"
        test_part = part.wald(formula_loc=formula,
                              factor_loc_totest="condition")
    else:
        formula = "~ 1 + condition + " + options.cov
        test_part = part.wald(formula_loc=formula,
                              factor_loc_totest="condition",
                              as_numeric=[options.cov])
    for i in clusters:
        df = test_part.tests[test_part.partitions.index(i)].summary()
        df.to_csv('{}-cluster_de_wald.tsv'.format(i), sep='\t', index=False)


if __name__ == '__main__':
    main()
