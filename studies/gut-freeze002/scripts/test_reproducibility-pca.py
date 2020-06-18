#!/usr/bin/env python

__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import os
import random
import numpy as np
import scanpy as sc
# import csv

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


def main():
    """Run CLI."""

    sc.settings.n_jobs = 30

    f1 = "/lustre/scratch119/realdata/mdt2/projects/sc-eqtl-ibd/analysis/leland_dev/nf_test/normalize=total_count.vars_to_regress=none.hvg_exclude=genes_remove_hvg_v001.scores=genes_score_v001/adata-normalized_pca.h5ad"
    f1 = "/lustre/scratch119/realdata/mdt2/projects/sc-eqtl-ibd/analysis/freeze_002/ti-cd_healthy/sample_qc=mito0pt80_ngene100_singlet-final_run_1pt8/nf-sample_qc=mito0pt80_ngene100_singlet-parameter_sweep_v002/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/adata-normalized_pca-bbknn.h5ad"

    # Load the AnnData file
    adata1 = sc.read_h5ad(filename=f1)
    adata2 = sc.read_h5ad(filename=f1)

    # Calculate pcs
    os.environ['PYTHONHASHSEED'] = str(seed_value)
    random.seed(seed_value)
    np.random.seed(seed_value)
    sc.tl.pca(
        adata1,
        n_comps=min(200, adata1.var['highly_variable'].sum()),
        zero_center=True,  # if false and arpack, not reproducible
        svd_solver='arpack',  # lobpcg not found in current sklearn
        use_highly_variable=True,
        copy=False,
        random_state=0,
        chunked=False
    )
    os.environ['PYTHONHASHSEED'] = str(seed_value)
    random.seed(seed_value)
    np.random.seed(seed_value)
    sc.tl.pca(
        adata2,
        n_comps=min(200, adata2.var['highly_variable'].sum()),
        zero_center=False,  # if false and arpack, not reproducible
        svd_solver='arpack',  # Scanpy default arpack as of 1.4.5
        use_highly_variable=True,
        copy=False,
        random_state=0,
        chunked=False
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
