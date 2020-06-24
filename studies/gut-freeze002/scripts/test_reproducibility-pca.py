#!/usr/bin/env python

__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import os
import random
import numpy as np
import scipy as sp
import scanpy as sc
import sklearn.utils
import sklearn.decomposition
# import csv

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


# def glm_pca(
#     data,
#     n_comps=None,
#     svd_solver='arpack',
#     use_highly_variable=None,
#     copy=False
# ):
#     """Runs GLM PCA
#     (1)
#
#     (2) Run GLM-PCA directly on UMI counts (or cp10k) using Poisson model
#     """
#     # Use Poisson model by default
#     glmpca()

def pca(
    data,
    n_comps=None,
    svd_solver='arpack',
    use_highly_variable=None,
    copy=False
):
    """Compute PCA coordinates, loadings and variance decomposition.

    Derived from scanpy 1.5.1.
    Principal component analysis [Pedregosa11]_.]
    Uses the implementation of *scikit-learn* [Pedregosa11]_.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    n_comps
        Number of principal components to compute. Defaults to 50, or 1 -
        minimum dimension size of selected representation.
    svd_solver
        SVD solver to use:
        `'arpack'` (the default)
          for the ARPACK wrapper in SciPy (:func:`~scipy.sparse.linalg.svds`)
        `'randomized'`
          for the randomized algorithm due to Halko (2009).
        `'auto'`
          chooses automatically depending on the size of the problem.
        `'lobpcg'`
          An alternative SciPy solver.
        .. versionchanged:: 1.4.5
           Default value changed from `'auto'` to `'arpack'`.
        Efficient computation of the principal components of a sparse matrix
        currently only works with the `'arpack`' or `'lobpcg'` solvers.
    use_highly_variable
        Whether to use highly variable genes only, stored in
        `.var['highly_variable']`.
        By default uses them if they have been determined beforehand.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.
    Returns
    -------
    adata : anndata.AnnData
        …otherwise if `copy=True` it returns or else adds fields to `adata`:
        `.obsm['X_pca']`
             PCA representation of data.
        `.varm['PCs']`
             The principal components containing the loadings.
        `.uns['pca']['variance_ratio']`
             Ratio of explained variance.
        `.uns['pca']['variance']`
             Explained variance, equivalent to the eigenvalues of the
             covariance matrix.
    """
    adata = data.copy() if copy else data

    if use_highly_variable and 'highly_variable' not in adata.var.keys():
        raise ValueError(
            'Did not find adata.var[\'highly_variable\']. '
            'Either your data already only consists of highly-variable genes '
            'or consider running `pp.highly_variable_genes` first.'
        )
    if use_highly_variable is None:
        if 'highly_variable' in adata.var.keys():
            use_highly_variable = True
        else:
            use_highly_variable = False

    if use_highly_variable:
        adata_comp = (
            adata[:, adata.var['highly_variable']]
        )
    else:
        adata_comp = adata

    if n_comps is None:
        min_dim = min(adata_comp.n_vars, adata_comp.n_obs)
        n_comps = min_dim - 1

    # random_state = sklearn.utils.check_random_state(random_state)
    X = adata_comp.X

    # If sparse, make dense.
    # Another option:
    # output = _pca_with_sparse(
    #     X, n_comps, solver=svd_solver, random_state=random_state
    # )
    if sp.sparse.issparse(X):
        X = X.toarray()

    # Sort out the solver
    if svd_solver == 'auto':
        svd_solver = 'arpack'
    if svd_solver not in {'arpack', 'randomized'}:
        raise ValueError(
            'svd_solver: {svd_solver} can not be used with sparse input.'
        )

    pca_ = sklearn.decomposition.PCA(
        n_components=n_comps,
        svd_solver=svd_solver,
        random_state=0
    )
    X_pca = pca_.fit_transform(X)

    # Cast to whatever datatype.
    # dtype = 'float32'
    # dtype
    #     Numpy data type string to which to convert the result.
    # if X_pca.dtype.descr != np.dtype(dtype).descr:
    #     X_pca = X_pca.astype(dtype)

    # Update the adata frame (if copy=False, then this is the same input adata
    # that the user provided)
    adata.obsm['X_pca'] = X_pca
    adata.uns['pca'] = {}
    adata.uns['pca']['params'] = {
        'zero_center': True,
        'use_highly_variable': use_highly_variable,
    }
    if use_highly_variable:
        adata.varm['PCs'] = np.zeros(shape=(adata.n_vars, n_comps))
        adata.varm['PCs'][adata.var['highly_variable']] = pca_.components_.T
    else:
        adata.varm['PCs'] = pca_.components_.T
    adata.uns['pca']['variance'] = pca_.explained_variance_
    adata.uns['pca']['variance_ratio'] = pca_.explained_variance_ratio_

    return adata if copy else None


def main():
    """Run CLI."""

    sc.settings.n_jobs = 30

    # Use pbmc3k dataset
    adata1 = sc.datasets.pbmc3k()
    sc.pp.filter_genes(adata1, min_counts=1)
    sc.pp.log1p(adata1)
    sc.pp.normalize_total(adata1)
    sc.pp.highly_variable_genes(adata1)

    f1 = "/lustre/scratch119/realdata/mdt2/projects/sc-eqtl-ibd/analysis/freeze_002/ti-cd_healthy/mito0pt80_ngene100_singlet-rep_test_09-multiplet_arpack/nf-sample_qc=mito0pt80_ngene100_singlet-parameter_sweep_v002/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/cluster.number_neighbors=-1.method=leiden.resolution=1.0/adata-normalized_pca-bbknn-umap-clustered.h5ad"
    adata1 = sc.read_h5ad(f1)
    adata2 = adata1.copy()

    # Calculate pcs
    # sc.tl.pca(
    #     adata1,
    #     n_comps=min(200, adata1.var['highly_variable'].sum()),
    #     zero_center=True,  # if false and arpack, not reproducible
    #     svd_solver='arpack',  # lobpcg not found in current sklearn
    #     use_highly_variable=True,
    #     copy=False,
    #     random_state=0,
    #     chunked=False
    # )
    pca(
        adata1,
        n_comps=min(200, adata1.var['highly_variable'].sum()),
        svd_solver='arpack',  # lobpcg not found in current sklearn
        use_highly_variable=True,
        copy=False
    )
    pca(
        adata2,
        n_comps=min(200, adata1.var['highly_variable'].sum()),
        svd_solver='arpack',  # lobpcg not found in current sklearn
        use_highly_variable=True,
        copy=False
    )

    dif_pc = np.isclose(adata1.obsm['X_pca'], adata2.obsm['X_pca'], atol=0)
    print((dif_pc == False).sum())
    print(adata1.obsm['X_pca'][np.invert(dif_pc)])
    print(adata2.obsm['X_pca'][np.invert(dif_pc)])






    # Add slight noise to adata2 PCS
    noise = np.random.normal(0, 1e-05, adata2.obsm['X_pca'].shape)
    adata2.obsm['X_pca'] = adata2.obsm['X_pca'] + noise

    # Run bbknn
    # Total number of neighbours = neighbors_within_batch x the number of
    # batches.
    os.environ['PYTHONHASHSEED'] = str(seed_value)
    random.seed(seed_value)
    np.random.seed(seed_value)
    sc.external.pp.bbknn(
        adata=adata1,
        batch_key='experiment_id',
        copy=False,
        # neighbors_within_batch=5
        n_pcs=29
    )
    os.environ['PYTHONHASHSEED'] = str(seed_value)
    random.seed(seed_value)
    np.random.seed(seed_value)
    sc.external.pp.bbknn(
        adata=adata2,
        batch_key='experiment_id',
        copy=False,
        # neighbors_within_batch=5
        n_pcs=29
    )

    # Compare the neighbor graph
    print(adata1.obsp['distances'].sum())
    print(adata2.obsp['distances'].sum())
    dist1 = adata1.obsp['distances'].todense()
    dist2 = adata2.obsp['distances'].todense()
    dif_dist = np.isclose(dist1, dist2)
    print((dif_dist == False).sum())
    print(dist1[np.invert(dif_dist)])
    print(dist2[np.invert(dif_dist)])

    print(adata1.obsp['connectivities'].sum())
    print(adata2.obsp['connectivities'].sum())
    conn1 = adata1.obsp['connectivities'].todense()
    conn2 = adata2.obsp['connectivities'].todense()
    dif_con = np.isclose(dist1, dist2)
    print((dif_con == False).sum())
    print(conn1[np.invert(dif_con)])
    print(conn2[np.invert(dif_con)])

    # UMAP
    # Saved to adata.uns['umap'] and adata.obsm['X_umap']
    sc.tl.umap(
        adata1,
        min_dist=0.05,  # Scanpy default = 0.05
        spread=1.0,  # Scanpy default = 1.0
        init_pos='spectral',  # Scanpy default = spectral
        n_components=29,
        # For some reason cannot access neighbors key slot, thus we
        # must keep uns['neighbors'] until we have run this.
        # neighbors_key='neighbors__{}'.format(plt__label),
        copy=False,
        random_state=0
    )
    sc.tl.umap(
        adata2,
        min_dist=0.05,  # Scanpy default = 0.05
        spread=1.0,  # Scanpy default = 1.0
        init_pos='spectral',  # Scanpy default = spectral
        n_components=29,
        # For some reason cannot access neighbors key slot, thus we
        # must keep uns['neighbors'] until we have run this.
        # neighbors_key='neighbors__{}'.format(plt__label),
        copy=False,
        random_state=0
    )

    # Compare UMAPS
    dif_umap = np.isclose(
        adata1.obsm['X_umap'],
        adata2.obsm['X_umap'],
        atol=1e-1
    )
    print((dif_umap == False).sum())
    print(adata1.obsm['X_umap'][np.invert(dif_umap)])
    print(adata2.obsm['X_umap'][np.invert(dif_umap)])


if __name__ == '__main__':
    main()
