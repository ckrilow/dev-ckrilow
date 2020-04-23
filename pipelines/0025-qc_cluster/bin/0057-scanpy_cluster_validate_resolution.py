#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
import scanpy as sc
import csv
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import metrics
from sklearn.linear_model import LogisticRegression
import joblib  # for numpy matrix, joblib faster than pickle


# based on NaiveDE
def create_colors(lr):
    n_cts = lr.classes_.shape[0]
    color_norm = colors.Normalize(vmin=-n_cts / 3, vmax=n_cts)
    ct_arr = np.arange(n_cts)
    ct_colors = cm.YlOrRd(color_norm(ct_arr))

    return ct_colors

# based on NaiveDE
def plot_roc(y_prob, y_test, lr):
    ct_colors = create_colors(lr)

    for i, cell_type in enumerate(lr.classes_):
        fpr, tpr, _ = metrics.roc_curve(y_test == cell_type, y_prob[:, i])
        plt.plot(fpr, tpr, c=ct_colors[i], lw=2)

    plt.plot([0, 1], [0, 1], color='k', ls=':')
    plt.xlabel('FPR')
    plt.ylabel('TPR')


# based on NaiveDE
def logistic_model(
    X,
    cell_labels,
    sparsity=1.0,
    test_size=0.25,
    n_jobs=None,
    verbose=True
):
    # Standarize features - this is especially important for SAGA
    # This scaler can also be applied to sparse CSR or CSC matrices by
    # passing with_mean=False to avoid breaking the sparsity structure
    # of the data.
    scaler = preprocessing.StandardScaler(
        with_mean=False,
        with_std=True
    )
    X_std = scaler.fit_transform(X)
    # X_std = (X / X.std()).dropna(1)  # Does not handle sparse matrix
    if verbose:
        print('Scaled X')
        print(type(X_std))

    # Split the data into training and test data.
    # NOTE: does it make sense to stratify sampling by cell_lables?
    X_train, X_test, y_train, y_test = model_selection.train_test_split(
        X_std,
        cell_labels,
        test_size=test_size
    )
    if verbose:
        print(
            'Split X into training {} and test {} sets.'.format(
                X_train.shape,
                X_test.shape
            )
        )
        print(type(X_train))

    # As noted in scanpy documentation:
    # penalty='l1' to try to come up with a minimal set of genes that are good
    # predictors (sparse solution meaning few non-zero fitted coefficients).
    #
    # SAGA is a good solver for large datasets and supports l1 penalty.
    # https://stackoverflow.com/questions/38640109/logistic-regression-python-solvers-defintions
    lr = LogisticRegression(
        penalty='l1',
        solver='saga',  # ‘liblinear’ and ‘saga’ also handle L1 penalty
        C=sparsity,
        n_jobs=n_jobs
    )
    lr.fit(X_train, y_train)
    if verbose:
        print('Completed: LogisticRegression.')

    y_prob = lr.predict_proba(X_test)
    if verbose:
        print('Completed: predict_proba.')
        print(type(y_prob))

    if verbose:
        print((lr.coef_ > 0).sum(1))

    lr_res = pd.DataFrame.from_records(
        lr.coef_
        # columns=X_std.columns
    )

    return y_prob, y_test, lr_res, lr


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and PCs file. Clusters the data. Saves an
            AnnData object with clusters in the clusters slot, a clusters
            file, and QC plots.
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
        help='H5 AnnData file where clusters have been saved to cluster slot.'
    )

    parser.add_argument(
        '-s', '--sparsity',
        action='store',
        dest='sparsity',
        default=0.25,
        type=float,
        help='Sparsity.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ts', '--test_size',
        action='store',
        dest='test_size',
        default=0.33,
        type=float,
        help='Test size fraction.\
            (default: %(default)s)'
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
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <h5_anndata>)'
    )

    options = parser.parse_args()

    # Fixed settings.
    # verbose = True

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}'.format(
            os.path.basename(options.h5.rstrip('.h5ad'))
        )

    # Load the AnnData file.
    # This file should already have clusters identified and saved to the
    # clusters slot.
    adata = sc.read_h5ad(filename=options.h5)

    # TODO: Optionally downsample cells

    # Set X to raw counts
    X = adata.layers['counts']
    # Alternative: set X to log1p_cp10k
    # X = adata.layers['log1p_cp10k']

    # NOTE: no need to make dense with sklearn... can keep spase
    # X = pd.DataFrame(X.toarray())

    # Here sparsity or C is the C param from
    # sklearn.linear_model.LogisticRegression
    # Inverse of regularization strength; must be a positive float.
    # Like in support vector machines, smaller values specify
    # stronger regularization.
    #
    # NOTE: for train / test split in logistic_model, this is done
    # via sklearn.model_selection.train_test_split... it may make
    # sense to set stratify to cell types.
    # cell_labels.logistic_model(
    test_prob, test_truth, lr_res, lr = logistic_model(
        X,
        adata.obs['cluster'].values,
        sparsity=options.sparsity,
        test_size=options.test_size,
        n_jobs=options.ncpu,
        verbose=True
    )

    # Save the lr_res dataframe.
    lr_res.to_csv(
        '{}-lr_res.tsv.gz'.format(out_file_base),
        sep='\t',
        index=True,
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression='gzip'
    )

    # Save the test results
    test_results = {
        'prob': test_prob,
        'truth': test_truth
    }
    joblib.dump(
        test_results,
        '{}-test_results.joblib.gz'.format(out_file_base),
        compress=('gzip', 5)
    )

    # Save the model
    joblib.dump(
        lr,
        '{}-lr.joblib.gz'.format(out_file_base),
        compress=('gzip', 5)
    )

    # Save the ROC of the test and truth.
    fig, ax = plt.subplot(2, 1, 2)
    plot_roc(test_prob, test_truth, lr)
    fig.savefig(
        '{}-roc.pdf'.format(out_file_base),
        dpi=300,
        bbox_inches='tight'
    )


if __name__ == '__main__':
    main()
