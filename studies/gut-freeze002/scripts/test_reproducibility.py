#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-05-29'
__version__ = '0.0.1'

import scanpy as sc
# import pandas as pd
import numpy as np
# import os


def run():
    """Run."""

    # Attempt 1
    #f1 = "/lustre/scratch119/realdata/mdt2/projects/sc-eqtl-ibd/analysis/freeze_002/ti-cd_healthy/sample_qc=mito0pt80_ngene100_singlet-final_run_1pt7/nf-sample_qc=mito0pt80_ngene100_singlet-parameter_sweep_v002/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/adata-normalized_pca-bbknn.h5ad"
    #f2 = "/lustre/scratch119/realdata/mdt2/projects/sc-eqtl-ibd/analysis/freeze_002/ti-cd_healthy/sample_qc=mito0pt80_ngene100_singlet-final_run_1pt8/nf-sample_qc=mito0pt80_ngene100_singlet-parameter_sweep_v002/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/adata-normalized_pca-bbknn.h5ad"
    # Attempt 2
    f1 = "/lustre/scratch119/realdata/mdt2/projects/sc-eqtl-ibd/analysis/freeze_002/ti-cd_healthy/mito0pt80_ngene100_singlet-rep_test_01/nf-sample_qc=mito0pt80_ngene100_singlet-parameter_sweep_v002/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/adata-normalized_pca-bbknn.h5ad"
    f2 = "/lustre/scratch119/realdata/mdt2/projects/sc-eqtl-ibd/analysis/freeze_002/ti-cd_healthy/mito0pt80_ngene100_singlet-rep_test_02/nf-sample_qc=mito0pt80_ngene100_singlet-parameter_sweep_v002/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/adata-normalized_pca-bbknn.h5ad"

    # load two runs
    adata1 = sc.read_h5ad(f1)
    adata2 = sc.read_h5ad(f2)

    # check that scrublet results are the same
    if np.array_equal(adata1.obs.index, adata2.obs.index) is False:
        print("filtering results different - check scrublet")

    # check that normalized expression is the same
    # dif_X = np.isclose(adata1.X, adata2.X, atol=0)
    # adata1.X[np.invert(dif_X)]
    # adata2.X[np.invert(dif_X)]

    # adata1.layers['log1p_cp10k']
    # adata2.layers['log1p_cp10k']

    # check that highly_variable genes are the same
    dif_hvg = np.equal(
        adata1.var['highly_variable'],
        adata2.var['highly_variable']
    )
    print((dif_hvg == False).sum())

    # check that PCs are the same
    X_pca1 = np.around(adata1.obsm['X_pca'], decimals=2)
    X_pca2 = np.around(adata2.obsm['X_pca'], decimals=2)
    dif_pc_v1 = np.equal(X_pca1, X_pca2)
    # The differences are exteremely minor - like 0.3 vs 0.4
    X_pca1[np.invert(dif_pc_v1)]
    X_pca2[np.invert(dif_pc_v1)]
    # The PCs are essentially the same:
    dif_pc = np.isclose(adata1.obsm['X_pca'], adata2.obsm['X_pca'], atol=1e-1)
    print((dif_pc == False).sum())
    print(adata1.obsm['X_pca'][np.invert(dif_pc)])
    print(adata2.obsm['X_pca'][np.invert(dif_pc)])
    if (dif_pc is False).sum() > 100:
        print("PCs may be different")

    # check that bbknn results are the same
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


if __name__ == '__main__':
    run()
