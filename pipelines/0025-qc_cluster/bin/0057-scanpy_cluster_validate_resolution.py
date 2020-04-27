#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-04-24'
__version__ = '0.0.1'

import argparse
import os
import numpy as np
import scipy as sci
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
import scanpy as sc
import csv
import joblib  # for numpy matrix, joblib faster than pickle

from sklearn import metrics
from sklearn import preprocessing
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression

# Import dask things that can be imported here.
from dask_jobqueue import LSFCluster
from dask.distributed import Client, performance_report, get_task_stream
# import dask
# import dask.array as da
# from dask_ml import preprocessing
# from dask_ml import model_selection
# import dask_ml.linear_model as dlm

# NOTE: The whole clustering pipeline could be replaced using dask as
# described here. One could do a grid search over (1) neighbors for clustering,
# (2) resolution for clustering, (3) regularization ... using the logreg score
# as the estimator to guide selection. See more below:
# https://ml.dask.org/examples/dask-glm.html


def start_dask_lsfcluster(cluster_size=5):
    """Start a dask cluster."""
    if cluster_size < 4:
        raise Exception('Too small of a cluster')
    # Settings for Sanger farm
    memory_in_gb = 20
    cluster = LSFCluster(
        queue='normal',
        walltime='00:30',
        log_directory='{}/dask_logs'.format(os.getcwd()),
        cores=4,
        memory='{} Gb'.format(memory_in_gb),
        mem=memory_in_gb*1e+9,  # should be in bytes
        lsf_units='mb',
        job_extra=[
            '-G team152',
            '-g /lt9/dask',
            '-R "select[mem>{}] rusage[mem={}]"'.format(
                int(memory_in_gb*1e+3),
                int(memory_in_gb*1e+3)
            )
        ],
        use_stdin=True
    )

    # View the job submission from Dask
    # cluster.job_script()

    # Scale cluster
    cluster.scale(cluster_size)

    # auto-scale between 10 and 100 jobs
    # cluster.adapt(
    #     minimum_jobs=int(cluster_size/4),
    #     maximum_jobs=cluster_size
    # )
    # cluster.adapt(maximum_memory="10 TB")  # use core/memory limits

    client = Client(
        cluster,
        timeout=120
    )
    client.wait_for_workers(n_workers=cluster_size)
    # print(client.scheduler_info()['services'])

    return cluster, client


def _create_colors(lr):
    n_cts = lr.classes_.shape[0]
    color_norm = colors.Normalize(vmin=-n_cts / 3, vmax=n_cts)
    ct_arr = np.arange(n_cts)
    ct_colors = cm.YlOrRd(color_norm(ct_arr))

    return ct_colors


def plot_roc(y_prob, y_test, lr):
    """Plot ROC curve. Based off of NaiveDE library."""
    ct_colors = _create_colors(lr)

    for i, cell_type in enumerate(lr.classes_):
        fpr, tpr, _ = metrics.roc_curve(y_test == cell_type, y_prob[:, i])
        plt.plot(fpr, tpr, c=ct_colors[i], lw=2)

    plt.plot([0, 1], [0, 1], color='k', ls=':')
    plt.xlabel('FPR')
    plt.ylabel('TPR')


def logistic_model(
    X,
    cell_labels,
    sparsity=1.0,
    test_size=0.25,
    n_jobs=None,
    standarize=True,
    with_mean=False,
    verbose=True,
    use_dask=False,
    out_file='out_file'
):
    """Fit logistic regression model. Based off of NaiveDE."""
    # Standarize features - this is especially important for SAGA
    # https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
    #
    # Many elements used in the objective function of a learning algorithm
    # (such as the RBF kernel of Support Vector Machines or the L1 and L2
    # regularizers of linear models) assume that all features are centered
    # around 0 and have variance in the same order. If a feature has a
    # variance that is orders of magnitude larger that others, it might
    # dominate the objective function and make the estimator unable to
    # learn from other features correctly as expected.
    #
    # NOTE: This scaler can also be applied to sparse CSR or CSC matrices by
    # passing with_mean=False to avoid breaking the sparsity structure
    # of the data.
    #
    # Standarization is also recommended by authors of
    # The Elements of Statistical Learning: Data Mining, Inference, and
    # Prediction. 9780387216065
    if standarize:
        scaler = preprocessing.StandardScaler(
            with_mean=with_mean,
            with_std=True
        )
        if with_mean and sci.sparse.issparse(X):
            X = X.todense()
        X_std = scaler.fit_transform(X)
        # X_std = (X / X.std()).dropna(1)  # Does not handle sparse matrix
        if verbose:
            print('Scaled X, with_mean={}'.format(with_mean))
    else:
        X_std = X

    # if use_dask:
    #     # Your chunks input will be normalized and stored in the third and
    #     # most explicit form. Note: chunks stands for "chunk shape" rather
    #     # than "number of chunks", so specifying chunks=1 means that you
    #     # will have many chunks, each with exactly one element.
    #     #
    #     # Automatic chunking
    #     # -1: no chunking along this dimension
    #     # None: no change to the chunking along this dimension (useful for
    #     #       rechunk)
    #     # "auto": allow the chunking in this dimension to accommodate ideal
    #     #         chunk sizes
    #     chunk_x = X.shape[0] // n_jobs
    #     chunk_y = X.shape[1] // n_jobs
    #     if verbose:
    #         print('Chunks: {},{}'.format(chunk_x, chunk_y))
    #     # NOTE: currently only allowed to chunk on x axis.
    #     X_std = da.from_array(X_std, chunks=({0: 'auto', 1: -1}))

    # Split the data into training and test data.
    # NOTE: does it make sense to stratify splitting by cell_lables?
    X_train, X_test, y_train, y_test = model_selection.train_test_split(
        X_std,
        cell_labels,
        # stratify=cell_labels,
        test_size=test_size
    )
    # if use_dask:
    #     # Calling dask.persist will preserve our data in memory, so no
    #     # computation will be needed as we pass over our data many times.
    #     # For example if our data came from CSV files and was not persisted,
    #     # then the CSV files would have to be re-read on each pass. This is
    #     # desirable if the data does not fit in RAM, but not slows down
    #     # our computation otherwise.
    #     X_train, X_test, y_train, y_test = dask.persist(
    #         X_train, X_test, y_train, y_test
    #     )
    if verbose:
        print(
            'Split X into training {} and test {} sets.'.format(
                X_train.shape,
                X_test.shape
            )
        )

    # As noted in scanpy documentation:
    # https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.rank_genes_groups.html
    # penalty='l1' to try to come up with a minimal set of genes that are good
    # predictors (sparse solution meaning few non-zero fitted coefficients).
    #
    # L1 = Lasso Regression, L2 = Ridge Regression. Because Lasso shrinks the
    # less important feature’s coefficient to zero (removing some feature
    # altogether), it works well for feature selection when there are a
    # large number of features.
    # https://scikit-learn.org/stable/auto_examples/linear_model/plot_logistic_l1_l2_sparsity.html
    #
    # SAGA is a good solver for large datasets  both the number of samples
    # and the number of features) and supports l1 penalty.
    # https://stackoverflow.com/questions/38640109/logistic-regression-python-solvers-defintions
    if use_dask:
        # lr = dlm.LogisticRegression(
        #     penalty='l1',
        #     C=sparsity
        # )
        lr = LogisticRegression(
            penalty='l1',
            solver='saga',  # ‘liblinear’ and ‘saga’ also handle L1 penalty
            C=sparsity,
            n_jobs=-1
        )
        # NOTE: Could add scater = [X_train] to joblib call to give a
        # local copy of X_train to each node.
        # NOTE: If no nodes are yet available in Client, code below will throw
        # error.
        with joblib.parallel_backend('dask'):
            with performance_report(
                filename='{}-dask-ml-performance_report.html'.format(
                    out_file
                )
            ):
                with get_task_stream(
                    filename='{}-dask-ml-task_stream.html'.format(
                        out_file
                    )
                ):
                    lr.fit(X_train, y_train)
    else:
        lr = LogisticRegression(
            penalty='l1',
            solver='saga',  # ‘liblinear’ and ‘saga’ also handle L1 penalty
            C=sparsity,
            n_jobs=n_jobs
        )
        lr.fit(X_train, y_train)
    if verbose:
        print('Completed: LogisticRegression.')
        print('Number of genes with > 0 coefficients:\t{}'.format(
            (lr.coef_ > 0).sum(1)
        ))

    # NOTE: P-value estimation
    # https://scikit-learn.org/stable/modules/linear_model.html#logistic-regression
    # It is possible to obtain the p-values and confidence intervals for
    # coefficients in cases of regression without penalization. The
    # statsmodels package <https://pypi.org/project/statsmodels/> natively
    # supports this. Within sklearn, one could use bootstrapping instead
    # as well.

    y_prob = lr.predict_proba(X_test)
    # if use_dask:
    #     y_prob = lr.predict_proba(X_test).compute()
    y_prob = pd.DataFrame(
        y_prob,
        columns=lr.classes_
    )
    y_prob['cell_label_true'] = y_test
    if verbose:
        print('Completed: predict_proba.')

    # Check out the performance of the model
    score_train = lr.score(X_train, y_train)
    score_test = lr.score(X_test, y_test)
    # if use_dask:
    #     score_train = score_train.compute()
    #     score_test = score_test.compute()
    if verbose:
        print('Training score:{}'.format(str(score_train)))
        print('Test score: {}'.format(str(score_test)))

    return lr, y_prob, score_train, score_test


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Fits logistic regression to predict labels in adata.obs["cluster"]'
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
        '-ncells', '--number_cells',
        action='store',
        dest='ncells',
        default=-1,
        type=int,
        help='Downsample to this number of cells.\
            (default: no downsampling)'
    )

    # NOTE: could potentially use a grid search estimator to set this
    # https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.GridSearchCV.html#sklearn.model_selection.GridSearchCV
    parser.add_argument(
        '-s', '--sparsity',
        action='store',
        dest='sparsity',
        default=0.1,
        type=float,
        help='LogisticRegression sparsity or inverse of regularization\
            strength; must be a positive float. Like in support vector\
            machines, smaller values specify stronger regularization.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ts', '--test_size',
        action='store',
        dest='test_size',
        default=0.33,
        type=float,
        help='Fraction of the data to use for test set.\
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
        '-d', '--dask_scale',
        action='store',
        dest='dask_scale',
        default=0,
        type=int,
        help='Scale using Dask if > 0.\
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
    verbose = True

    # Scanpy settings
    # sc.settings.figdir = os.getcwd()  # figure output directory to match base
    # sc.settings.n_jobs = options.ncpu  # number CPUs
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

    # Optionally downsample cells
    if options.ncells > -1:
        n_cells_start = adata.n_obs
        if n_cells_start < options.ncells:
            raise Exception('Fewer cells than ncells specified.')
        # Downsample cells by fraction
        # sc.pp.subsample(
        #     adata,
        #     fraction=0.2,
        #     copy=False
        # )
        # Downasmple cells to a specific number.
        sc.pp.subsample(
            adata,
            n_obs=options.ncells,
            copy=False
        )
        print('Cell downsample applied: {} dropped, {} remain.'.format(
            n_cells_start - adata.n_obs,
            adata.n_obs
        ))

    # Set X to log1p_cp10k
    adata.X = adata.layers['log1p_cp10k']
    # Set X to raw counts
    # adata.X = adata.layers['counts']

    # Ensure un-informative genes are filtered
    # sc.pp.filter_genes(adata, min_cells=5)

    # Subset to only highly variable genes
    subset_to_hcg = False
    if subset_to_hcg:
        adata = adata[:, adata.var['highly_variable']]

    # NOTE: We could use scanpy to stanardize the data, but just like sklearn
    #       this will also change the matrix from Compressed Sparse Row format
    #       to dense format.
    # sc.pp.scale(
    #     adata,
    #     zero_center=True,
    #     max_value=None,
    #     copy=False
    # )

    # NOTE: no need to make dense with sklearn... can keep spase
    # X = pd.DataFrame(X.toarray())

    use_dask = False
    n_jobs = options.ncpu
    if options.dask_scale > 0:
        cluster, client = start_dask_lsfcluster(options.dask_scale)
        use_dask = True
        n_jobs = int(options.dask_scale*1.25)
        print(client)

    # Here sparsity or C is the C param from
    # sklearn.linear_model.LogisticRegression
    # Inverse of regularization strength; must be a positive float.
    # Like in support vector machines, smaller values specify
    # stronger regularization.
    lr, test_results, score_train, score_test = logistic_model(
        adata.X,
        adata.obs['cluster'].values,
        sparsity=options.sparsity,
        test_size=options.test_size,
        n_jobs=n_jobs,
        standarize=True,
        with_mean=True,
        verbose=verbose,
        use_dask=use_dask,
        out_file=out_file_base
    )
    if use_dask:
        cluster.close()

    # Save the model
    out_f = '{}-lr_model.joblib.gz'.format(out_file_base)
    joblib.dump(
        lr,
        out_f,
        compress=('gzip', 3)
    )
    # Example of how to load a model.
    # lr = joblib.load(
    #     'CD5677E01F0E30D5-adata-normalized_pca-clustered-lr.joblib.gz'
    # )
    if verbose:
        print('Completed: save {}.'.format(out_f))

    # Save the test results - each row is a cell and the columns are the prob
    # of that cell belonging to a particular class.
    out_f = '{}-test_result.tsv.gz'.format(out_file_base)
    test_results.to_csv(
        out_f,
        sep='\t',
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression='gzip'
    )
    if verbose:
        print('Completed: save {}.'.format(out_f))

    # Save the ROC of the test and truth.
    out_f = '{}-roc.pdf'.format(out_file_base)
    fig = plt.figure()
    cell_label_true = test_results.pop('cell_label_true')
    plot_roc(test_results.values, cell_label_true.values, lr)
    fig.savefig(
        out_f,
        dpi=300,
        bbox_inches='tight'
    )
    if verbose:
        print('Completed: save {}.'.format(out_f))

    # Make a dataframe of the coefficients.
    # Rows = cell type label and columns = genes tested
    lr_res = pd.DataFrame.from_records(
        lr.coef_,
        index=lr.classes_,
        columns=adata.var.index
    )
    # Save the lr_res dataframe.
    out_f = '{}-lr_coef.tsv.gz'.format(out_file_base)
    lr_res.to_csv(
        out_f,
        sep='\t',
        index=True,
        index_label='cell_label',
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression='gzip'
    )
    if verbose:
        print('Completed: save {}.'.format(out_f))


if __name__ == '__main__':
    main()
