#!/path/to/celltypist/python
# -*- coding: UTF-8 -*-
# # # # # # # # # # # # 
"""
╭───────────────────────────────────────╮ 
│ classify.py   2025/01/16/-16:42 
╰───────────────────────────────────────╯ 
│ Description:
    CellTypist Classify with trained model
""" # [By: Yiwen Hu & Synni Meng]

# region |- Tissue Only Celltype Dict -|
LIMIT_TISSUE = {
'adipose tissue': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'adipocytes', 'dendritic cells', 'smooth muscle cells', 'fibroblasts', 'monocytes', 'endothelial cells', 'schwann cells', 'b-cells'],
'bone marrow': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'erythroid cells', 'monocytes', 'dendritic cells', 'b-cells'],
'brain': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'oligodendrocyte precursor cells', 'microglial cells', 'fibroblasts', 'monocytes', 'dendritic cells', 'oligodendrocytes', 'endothelial cells', 'inhibitory neurons', 'excitatory neurons', 'b-cells', 'astrocytes'],
'breast': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'breast myoepithelial cells', 'adipocytes', 'dendritic cells', 'smooth muscle cells', 'fibroblasts', 'breast glandular cells', 'monocytes', 'endothelial cells', 'b-cells'],
'bronchus': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'basal respiratory cells', 'ionocytes', 'fibroblasts', 'monocytes', 'dendritic cells', 'club cells', 'endothelial cells', 'smooth muscle cells', 'b-cells', 'ciliated cells'],
'colon': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'enteroendocrine cells', 'distal enterocytes', 'fibroblasts', 'monocytes', 'dendritic cells', 'intestinal goblet cells', 'endothelial cells', 'smooth muscle cells', 'undifferentiated cells', 'b-cells'],
'endometrium': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'endometrial stromal cells', 'fibroblasts', 'monocytes', 'dendritic cells', 'glandular and luminal cells', 'endothelial cells', 'smooth muscle cells', 'b-cells', 'ciliated cells'],
'esophagus': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'basal squamous epithelial cells', 'dendritic cells', 'smooth muscle cells', 'fibroblasts', 'squamous epithelial cells', 'monocytes', 'endothelial cells', 'schwann cells', 'b-cells'],
'eye': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'cone photoreceptor cells', 'fibroblasts', 'rod photoreceptor cells', 'monocytes', 'dendritic cells', 'horizontal cells', 'endothelial cells', 'smooth muscle cells', 'b-cells', 'bipolar cells', 'muller glia cells'],
'fallopian tube': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'neutrophils', 'dendritic cells', 'secretory cells', 'smooth muscle cells', 'lymphatic endothelial cells', 'fibroblasts', 'monocytes', 'endothelial cells', 'b-cells', 'ciliated cells'],
'heart muscle': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'fibroblasts', 'monocytes', 'dendritic cells', 'cardiomyocytes', 'endothelial cells', 'smooth muscle cells', 'b-cells'],
'kidney': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'collecting duct cells', 'fibroblasts', 'distal tubular cells', 'monocytes', 'dendritic cells', 'endothelial cells', 'smooth muscle cells', 'proximal tubular cells', 'b-cells'],
'liver': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'erythroid cells', 'dendritic cells', 'smooth muscle cells', 'kupffer cells', 'hepatocytes', 'fibroblasts', 'monocytes', 'endothelial cells', 'b-cells', 'cholangiocytes'],
'lung': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'dendritic cells', 'alveolar cells type 1', 'smooth muscle cells', 'fibroblasts', 'monocytes', 'club cells', 'alveolar cells type 2', 'endothelial cells', 'b-cells', 'ciliated cells'],
'lymph node': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'fibroblasts', 'monocytes', 'dendritic cells', 'endothelial cells', 'smooth muscle cells', 'b-cells'],
'ovary': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'fibroblasts', 'monocytes', 'dendritic cells', 'granulosa cells', 'ovarian stromal cells', 'endothelial cells', 'smooth muscle cells', 'b-cells', 'lymphatic endothelial cells', 'oocytes'],
'pancreas': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'dendritic cells', 'exocrine glandular cells', 'smooth muscle cells', 'ductal cells', 'fibroblasts', 'monocytes', 'endothelial cells', 'pancreatic endocrine cells', 'b-cells'],
'pbmc': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'platelets', 'monocytes', 'dendritic cells', 'b-cells'],
'placenta': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'extravillous trophoblasts', 'syncytiotrophoblasts', 'dendritic cells', 'hofbauer cells', 'smooth muscle cells', 'cytotrophoblasts', 'fibroblasts', 'monocytes', 'endothelial cells', 'b-cells'],
'prostate': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'dendritic cells', 'smooth muscle cells', 'prostatic glandular cells', 'fibroblasts', 'monocytes', 'basal prostatic cells', 'club cells', 'endothelial cells', 'b-cells'],
'rectum': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'enteroendocrine cells', 'distal enterocytes', 'fibroblasts', 'monocytes', 'dendritic cells', 'intestinal goblet cells', 'endothelial cells', 'smooth muscle cells', 'paneth cells', 'undifferentiated cells', 'b-cells'],
'salivary gland': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'dendritic cells', 'smooth muscle cells', 'salivary duct cells', 'fibroblasts', 'mucus glandular cells', 'serous glandular cells', 'monocytes', 'endothelial cells', 'b-cells'],
'skeletal muscle': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'skeletal myocytes', 'fibroblasts', 'monocytes', 'dendritic cells', 'endothelial cells', 'smooth muscle cells', 'b-cells'],
'skin': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'suprabasal keratinocytes', 'langerhans cells', 'basal keratinocytes', 'dendritic cells', 'smooth muscle cells', 'fibroblasts', 'monocytes', 'endothelial cells', 'melanocytes', 'b-cells'],
'small intestine': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'enteroendocrine cells', 'fibroblasts', 'monocytes', 'dendritic cells', 'intestinal goblet cells', 'proximal enterocytes', 'endothelial cells', 'smooth muscle cells', 'paneth cells', 'undifferentiated cells', 'b-cells'],
'spleen': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'fibroblasts', 'monocytes', 'dendritic cells', 'endothelial cells', 'smooth muscle cells', 'b-cells'],
'stomach': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'fibroblasts', 'gastric mucus-secreting cells', 'monocytes', 'dendritic cells', 'endothelial cells', 'smooth muscle cells', 'b-cells'],
'testis': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'spermatogonia', 'dendritic cells', 'leydig cells', 'smooth muscle cells', 'spermatocytes', 'fibroblasts', 'late spermatids', 'monocytes', 'peritubular cells', 'endothelial cells', 'b-cells', 'early spermatids', 'sertoli cells'],
'thymus': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'dendritic cells', 'smooth muscle cells', 'fibroblasts', 'mesothelial cells', 'monocytes', 'endothelial cells', 'b-cells'],
'tongue': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'suprabasal keratinocytes', 'basal keratinocytes', 'dendritic cells', 'smooth muscle cells', 'fibroblasts', 'serous glandular cells', 'monocytes', 'endothelial cells', 'b-cells'],
'vascular': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'fibroblasts', 'monocytes', 'dendritic cells', 'endothelial cells', 'smooth muscle cells', 'schwann cells', 'b-cells'],
'glial cells': ["astrocytes", "oligodendrocyte precursor cells", "oligodendrocytes","microglial cells", "macrophages", "monocytes", "dendritic cells", "granulocytes", "t-cells", "nk-cells"],
'thyroid gland': ['plasma cells', 'nk-cells', 'granulocytes', 't-cells', 'macrophages', 'thyroid glandular cells', 'b-cells', 'dendritic cells', 'endothelial cells', 'fibroblasts', 'monocytes', 'smooth muscle cells'],
}
# endregion

# region |- Import -|
# from distutils import LooseVersion
from contextlib import contextmanager  
import celltypist
from celltypist import models
import scanpy as sc
import argparse
import time
import numpy as np
import pandas as pd
import os

import sklearn.metrics as metrics
# endregion

parser = argparse.ArgumentParser(description='CellTypist Classify with trained model')
parser.add_argument('-i', '--input', type=str, help='Input unClassified h5ad path')
parser.add_argument('-m', '--model', type=str, help='Trained Model path')
# Allow multiple tissues as input, separated by commas
parser.add_argument('-t', '--tissue', type=str, nargs="+", default=None, help='Tissue Limitation', required=False)
# parser.add_argument('-t', '--tissue', type=str, nargs='+', default=None, help='Tissue Limitation', required=False)

parser.add_argument('-o', '--output', type=str, help='Output Prefix')
parser.add_argument('-l', '--label', type=str, default=None, help='Label key in adata.obs.keys()')
parser.add_argument('-d', '--detailed', action="store_true", help='make Plots')
parser.add_argument('-nm', '--nmarkers', type=int, default=3, help='Plot top n markers')
args = parser.parse_args()


args.model = args.model if args.model.endswith(".pkl") else f"{args.model}/model.pkl"
print("Args:\n",
    f"Use Model: {args.model}\n",
    f"Tissue: [ {args.tissue} ]\n",
    f"Predict On: {args.input}\n",
    f"Output: {args.output}\n")

# region |- Functions -|
@contextmanager
def Duration(name='Job-Name', unit='s'):
    """
    Context manager that prints the execution time of the corresponding code block
    Args:
        name: Custom task name
        unit: Time unit—h (hours), m (minutes), or s (seconds)
    Returns:

    """
    unit = unit.lower()[0] 
    private_map = {'s': 1, 'm': 60, 'h': 3600}
    factor = private_map[unit]
    t = time.time()  # record time
    yield t     # time
    print(f"\n[ Elapsed time ] {name}：{(time.time()-t)/factor:.2f} {unit}")

# endregion

adata = sc.read_h5ad(args.input)
adata.var_names = adata.var_names.astype(str)

# is data normalized ?
if adata.X.max().is_integer():
    print("Preprocess data...")
    # norm
    sc.pp.normalize_total(adata, target_sum = 1e4)
    sc.pp.log1p(adata)
else:
    # check norm val
    exp_sum = adata.X.expm1().sum(axis = 1)
    if abs(exp_sum.mean() - 1e4) > 1 or exp_sum.var() > 0.1:
        raise ValueError('Input h5ad Counts is not normed with 1e4! Please use raw counts!')


model = models.Model.load(model = args.model)
# if tissue(s) are provided, use the corresponding cell types
# allow_cates = LIMIT_TISSUE[args.tissue] if args.tissue else model.cell_types
if args.tissue != None:
    # Validate tissue input, accounting for the shell replacing spaces with underscores
    allow_cates = []
    for tissue in args.tissue:
        tissue = tissue.lower().replace("_", " ")
        assert tissue in LIMIT_TISSUE.keys(), f"{tissue} is not in Tissue keys\n{LIMIT_TISSUE.keys()}"
        allow_cates.extend(LIMIT_TISSUE.get(tissue, []))  
    allow_cates = np.unique(allow_cates).tolist() 
else:
    allow_cates = model.cell_types  # If no tissue is provided, use all cell types available in the model

# Ensure the specified cell types are within the model’s predictable set
out_cates = set(allow_cates) - set(model.cell_types) 
assert len(out_cates) == 0, f"Error: {out_cates} is not in model prediction celltypes set!"


with Duration('Predict', unit='m'):
    predictions = celltypist.annotate(adata, model = model, majority_voting = False)
    pred_scores = predictions.probability_matrix    # 	Confidence scores range 0–1; per cell, probabilities sum to 1 across cell types
    # Assign the highest-probability cell type among the allowed ones
    adata.obs['predicted_labels'] = pred_scores[allow_cates].idxmax(axis=1) # 在允许被标注出的 细胞 上获得 最大概率的celltype
    # Assign the highest-probability cell type among the allowed ones
    max_scores = pred_scores[allow_cates].max(axis=1)
    # Save predicted labels and their scores to CSV
    result_df = pd.DataFrame({
        'predicted_labels': adata.obs['predicted_labels'],  # Predicted labels
        'scores': max_scores  # Maximum scores
    }, index=adata.obs.index)  # Set index to cell names
    result_df.to_csv(f"{args.output}.pred_labels_scores.csv", index=True)
    # # Conditionally assign: replace with "unassigned" if max_scores < n_max
    # n_max = 0.5
    # adata.obs['predicted_labels'] = adata.obs['predicted_labels'].where(max_scores >= n_max, "unassigned")
    # Save limited cell-type scores to CSV; float_format='%.6f' keeps 6 decimal places
    pred_scores[allow_cates].to_csv(f"{args.output}.pred_limit_scores.csv", float_format='%.6f')
    # Uncomment below if needed
    # pred_scores.to_csv(f"{args.output}.pred_scores.csv", float_format='%.6f')
    # Attach prediction scores to adata.uns (comment out the next two lines if not needed)
    adata.uns['pred_scores'] = pred_scores      # Full cell-type scores
    adata.uns['pred_limit_scores'] = adata.uns['pred_scores'][allow_cates]   # Scores restricted to allowed cell types
    adata.obs["predicted_labels"].to_csv(f"{args.output}.pred_labels.csv")


if args.detailed:
    if 'neighbors' not in adata.obsm: 
        print("compute hvg & neighbors...")
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
    print("Plot UMAP...")
    # Plot UMAP and heatmap
    if 'X_umap' not in adata.obsm: sc.tl.umap(adata)
    colors = ('predicted_labels', args.label) if args.label else 'predicted_labels'
    fig_umap = sc.pl.umap(adata, color = colors, return_fig=True, wspace=0.5)    # legend_loc = 'on data', # wspace
    fig_umap.savefig(f"{args.output}.pred_umap.png", bbox_inches='tight', dpi=300)
    
    # 	Add the max value of pred_limit_scores to adata.obs
    adata.obs['pred_limit_scores'] = adata.uns['pred_limit_scores'].max(axis=1)
    fig_score = sc.pl.umap(adata, color = "pred_limit_scores", color_map='viridis', return_fig=True, wspace=0.5)
    fig_score.savefig(f"{args.output}.pred_score_umap.png", bbox_inches='tight', dpi=300)


if args.label:
    print("Plot Label Transfer HeatMap...")
    fig_hm = celltypist.dotplot(predictions, use_as_reference = args.label, use_as_prediction = 'predicted_labels', return_fig=True)
    fig_hm.savefig(f"{args.output}.transfer.png", bbox_inches='tight', dpi=300)
    # metrics
    targ = adata.obs[args.label].values
    pred = adata.obs['predicted_labels'].values
    acc = round(100*metrics.accuracy_score(targ, pred), 2),    # acc
    f1 = round(100*metrics.f1_score(targ, pred, average='macro'), 2),  # f1 macro
    pd.DataFrame({'Accuracy': acc, 'F1-Macro': f1}, index=[0]).to_csv(f"{args.output}.metrics.csv", index=False, sep='\t')



print("Predict adata Success.")
