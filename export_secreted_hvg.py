import pandas as pd
from pathlib import Path

hvg_path = Path("HVG.csv")
secreted_path = Path("secreted_gene_names.csv")
out_path = Path("secreted_HVG.csv")

hvg = pd.read_csv(hvg_path)
if "gene_name" not in hvg.columns:
    raise KeyError(f"{hvg_path} 缺少 gene_name 列，实际列: {list(hvg.columns)}")

secreted_df = pd.read_csv(secreted_path, header=None, names=["gene_symbol"])
secreted_genes = secreted_df["gene_symbol"].astype(str).str.strip()
if not secreted_genes.empty and secreted_genes.iloc[0].lower() == "gene_symbol":
    secreted_genes = secreted_genes.iloc[1:]
secreted_genes = secreted_genes[secreted_genes != ""].unique()

mask = hvg["gene_name"].isin(secreted_genes)
secreted_hvg = hvg[mask].copy()

secreted_hvg.to_csv(out_path, index=False)
print(f"筛到 {mask.sum()} 个 secreted 基因的 HVG，已写入 {out_path}")



