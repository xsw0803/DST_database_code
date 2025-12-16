import os
import re
import anndata as ad
import scanpy as sc
import pandas as pd

os.environ.setdefault("NUMBA_CACHE_DIR", "./numba_cache")

data_path = "normalized_expression.h5ad"
adata = sc.read(data_path)
adata.var_names_make_unique()

gene_dict = pd.read_csv(
    "master_dictionary_gene_unique.csv",
    encoding="utf-8",
    encoding_errors="replace",
    engine="python",
    on_bad_lines="skip",
)
gene_list = (
    gene_dict["gene_symbol"]
    .dropna()
    .astype(str)
    .str.strip()
    .unique()
)
genes_in_data = [g for g in gene_list if g in adata.var_names]
if not genes_in_data:
    raise ValueError("master_dictionary_gene_unique.csv is mismatched with gene_symbol in master dictionary.")

adata = adata[:, genes_in_data].copy()
print(f"Match {len(genes_in_data)} genes for differential gene analysis.")

meta_candidates = [
    "GSE174367_snRNA-seq_cell_meta.csv",
    "PFC_metadata_whole.csv",
]
meta_path = next((p for p in meta_candidates if os.path.exists(p)), None)
if meta_path is None:
    raise FileNotFoundError("No meta_csv, please provide GSE174367_snRNA-seq_cell_meta.csv or PFC_metadata_whole.csv")

with open(meta_path, "rb") as fh:
    header = fh.read(4)
try:
    if header.startswith(b"PK") or meta_path.lower().endswith((".xlsx", ".xls")):
        meta = pd.read_excel(meta_path)
    else:
        meta = pd.read_csv(
            meta_path,
            encoding="utf-8",
            encoding_errors="replace",
            engine="python",
            on_bad_lines="skip",
        )
except Exception:
    meta = pd.read_excel(meta_path)


if "Barcode" not in meta.columns:
    raise KeyError(f"{meta_path} 中缺少 Barcode 列，当前列: {list(meta.columns)}")

meta = meta.set_index("Barcode")
adata = adata[adata.obs_names.isin(meta.index)].copy()
meta_aligned = meta.loc[adata.obs_names]
meta_to_add = meta_aligned[[c for c in meta_aligned.columns if c not in adata.obs.columns]]
adata.obs = adata.obs.join(meta_to_add, how="left")

if "Cell_type" not in adata.obs.columns:
    raise ValueError("注释中缺少 Cell_type 列，无法按细胞类型分组做差异分析。")

cell_types = sorted(adata.obs["Cell_type"].dropna().unique())

for ct in cell_types:
    ad_ct = adata[adata.obs["Cell_type"] == ct].copy()
    if ad_ct.n_obs < 20:
        print(f"Skip {ct}: Too less cell numbers ({ad_ct.n_obs})")
        continue

    sc.tl.rank_genes_groups(
        ad_ct,
        groupby="Diagnosis",
        reference="CT",
        method="wilcoxon",
    )
    sc.tl.filter_rank_genes_groups(
        ad_ct,
        min_in_group_fraction=0.1,
        max_out_group_fraction=0.1,
        min_fold_change=1.5,
    )

    deg_df = sc.get.rank_genes_groups_df(ad_ct, group=None)

    safe_ct = re.sub(r"[^0-9A-Za-z._-]+", "_", ct)
    out_path = f"DEG_{safe_ct}.csv"
    deg_df.to_csv(out_path, index=False)
    deg_h5ad = ad.AnnData(
        X=None,
        obs=deg_df.reset_index(drop=True),
        var=pd.DataFrame(index=[]),
        uns={
            "groupby": "Diagnosis",
            "reference": "CT",
            "method": "wilcoxon",
            "cell_type": ct,
            "n_cells": ad_ct.n_obs,
            "filter_params": {
                "min_in_group_fraction": 0.1,
                "max_out_group_fraction": 0.1,
                "min_fold_change": 1.5,
            },
        },
    )
    deg_h5ad.write(f"DEG_{safe_ct}.h5ad", compression="gzip")
    print(f"{ct}: {ad_ct.n_obs} cells, DEG wirite in {out_path}.")
