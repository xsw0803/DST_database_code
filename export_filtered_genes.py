import os
import pandas as pd

os.environ.setdefault("NUMBA_CACHE_DIR", "./numba_cache")
data_path = "normalized_expression.h5ad" if os.path.exists("normalized_expression.h5ad") else "Norm_HVG.h5ad"

import scanpy as sc

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
    raise ValueError("master_dictionary_gene_unique.csv is dismtached with gene_symbol in master dictionary.")

out_path = "secreted_gene_names.csv"
pd.DataFrame({"gene_symbol": genes_in_data}).to_csv(out_path, index=False)

print(f"Write {len(genes_in_data)} genes into {out_path}.")


