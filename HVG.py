import os
import scanpy as sc
import pandas as pd
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "./numba_cache")

if os.path.exists("normalized_expression.h5ad"):
    adata = sc.read("normalized_expression.h5ad")
else:
    h5_path = "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5"
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.write("normalized_expression.h5ad", compression="gzip")

adata.var_names_make_unique()

gene_path = Path("master_dictionary_gene_unique.csv")
if gene_path.exists():
    gene_dict = pd.read_csv(
        gene_path,
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
    if genes_in_data:
        adata = adata[:, genes_in_data].copy()
        print(f"过滤到 {len(genes_in_data)} 个 master_dictionary_gene_unique 基因用于 HVG。")
    else:
        raise ValueError("master_dictionary_gene_unique.csv 中的 gene_symbol 与表达矩阵不匹配。")

batch_key = None
if os.path.exists("GSE174367_snRNA-seq_cell_meta.csv"):
    meta_path = "GSE174367_snRNA-seq_cell_meta.csv"
    meta = None
    last_err = None
    for enc in [
        {"encoding": "utf-8", "encoding_errors": "replace"},
        {"encoding": "latin1"},
    ]:
        try:
            meta = pd.read_csv(
                meta_path,
                engine="python",
                on_bad_lines="skip",
                **enc,
            )
            break
        except UnicodeDecodeError as e:
            last_err = e
    if meta is None or "Barcode" not in meta.columns:
        try:
            meta = pd.read_excel(meta_path)
        except Exception as e:
            raise KeyError(
                f"{meta_path} 缺少 Barcode 列，实际列: {list(meta.columns) if meta is not None else '未读出'}。"
                f"CSV 解码错误: {last_err}; Excel 读取错误: {e}"
            )
    if "Barcode" not in meta.columns:
        raise KeyError(f"{meta_path} 缺少 Barcode 列，实际列: {list(meta.columns)}")
    meta = meta.set_index("Barcode")
    adata = adata[adata.obs_names.isin(meta.index)].copy()
    adata.obs = adata.obs.join(meta, how="left")
    if "Batch" in adata.obs.columns:
        batch_key = "Batch"

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat_v3",
    batch_key=batch_key,
)

hvg_df_raw = adata.var[adata.var["highly_variable"]].copy()
hvg_df = hvg_df_raw.assign(
    gene_name=hvg_df_raw.index,
    mean=hvg_df_raw["means"],
    variance=hvg_df_raw["variances"],
    variance_norm=hvg_df_raw["variances_norm"],
)[["gene_name", "mean", "variance", "variance_norm"]]
hvg_df.to_csv("HVG.csv", index=False)

adata_hvg = adata[:, adata.var["highly_variable"]].copy()
adata_hvg.write("Norm_HVG.h5ad", compression="gzip")

print(f"找到 {hvg_df.shape[0]} 个高变基因，已保存为HVG.csv 和 Norm_HVG.h5ad")
