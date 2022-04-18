import os.path
import numpy as np
import pandas as pd
import scipy
import scipy.io
from scipy.sparse import csr_matrix, csc_matrix
import anndata
import scanpy as sc
import sys

#os.chdir('~/data/10_SingleCell/PaperReproduce/Skin_disease/HCA_skin-MHskin1/haniffalab-HCA_skin-2f859da')
adata = sc.read_h5ad('submission.h5ad')



sc.pp.filter_cells(adata, min_genes=50)
sc.pp.filter_genes(adata, min_cells=3)
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1


#sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True)

adata = adata[adata.obs['n_genes'] < 6000, :]
adata = adata[adata.obs['n_genes'] > 400, :]
adata = adata[adata.obs['n_counts'] > 1000, :]
adata = adata[adata.obs['percent_mito'] < 0.2, :]

#adata # n_obs × n_vars = 473325 × 28680
#adata.raw = adata #For backup


sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, adata.X, min_mean=0.0125, max_mean=5, min_disp=0.25)
sc.pl.highly_variable_genes(adata)


filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.0125, max_mean=5, min_disp=0.25)
adata_filtered = adata[:, filter_result.gene_subset]

sc.tl.pca(adata_filtered, n_comps=50)

import bbknn
bbknn.bbknn(adata_filtered, batch_key='donor_id', copy=False)
save_file = 'adata_filtered.h5ad'
adata_filtered.write(save_file)

