#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# # # # # # # # # # # # 
"""
╭───────────────────────────────────────╮ 
│ ScanpyFilter.py   2025/03/04/-10:07 
╰───────────────────────────────────────╯ 
│ Description:
    
""" # [By: HuYw]

# region |- Import -|
import scanpy as sc
import scipy as sp
import pandas as pd
import argparse
# endregion



parser = argparse.ArgumentParser(description="Process STAR Solo output for Scanpy Filter analysis.")
parser.add_argument("-i", "--input", required=True, type=str, help="Path to the STAR Solo output directory.")
parser.add_argument("-o", "--outdir", required=True, type=str, help="Output dir for the processed AnnData object (h5ad format).")
parser.add_argument("-j", "--job", required=True, type=str, help="Job name for output files.")
parser.add_argument("-c", "--core", required=False, default=4, type=int, help="n CPU core can use.")
args = parser.parse_args()
# n job / core
sc.settings.n_jobs = args.core

def readSTARsolo(path="."):
    barcodes = pd.read_csv(f"{path}/barcodes.tsv", sep="\t", index_col=0, header=None)
    features = pd.read_csv(f"{path}/features.tsv", sep="\t", index_col=0, header=None)
    features.columns = ("Symbol", "Feature Type")
    features['Feature ID'] = features.index.values
    features.index = features["Symbol"].copy()
    mtx = sp.io.mmread(f"{path}/matrix.mtx").T
    return sc.AnnData(X=sp.sparse.csr_matrix(mtx), 
                      obs=barcodes, var=features)

# 读取 STAR solo out
adata = readSTARsolo(args.input)
adata.var_names_make_unique()
# scanpy output Summary UMI: total_counts; nGene: n_genes_by_counts;
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], inplace=True, log1p=False
)


# filter cells
sc.pp.scrublet(adata)
adata = adata[
    (adata.obs['n_genes_by_counts'] > 500) &
    (adata.obs['total_counts'] > 1000) &
    (adata.obs['pct_counts_mt'] < 10) &
    ~adata.obs['predicted_doublet']
]

# preprocess
adata.layers["counts"] = adata.X.copy() # Saving count data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

# output
adata.write(f"{args.outdir}/{args.job}.h5ad")
pd.DataFrame(adata.obs_names).to_csv(f"{args.outdir}/{args.job}.cells.csv", header=False, index=False)

