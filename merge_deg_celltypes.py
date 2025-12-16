import anndata as ad
import pandas as pd
from pathlib import Path

deg_files = [
    fp
    for fp in sorted(Path(".").glob("DEG_*.csv"))
    if fp.stem not in {"DEG_celltype_merged", "DEG_Log_CellType"}
]
if not deg_files:
    raise FileNotFoundError("NO DEG_*.csv file，please run DEG_celltype.py")

merged = []
for fp in deg_files:
    try:
        df = pd.read_csv(fp)
    except Exception as exc:
        raise RuntimeError(f"Reading {fp} faili: {exc}") from exc
    
    drop_cols = {"Gene", "ASC", "EX", "INH", "MG", "ODC", "CellType", "OPC", "PER.END", "cell_type"}
    df = df.drop(columns=[c for c in drop_cols if c in df.columns], errors="ignore")

    cell_type = fp.stem.removeprefix("DEG_")
    df.insert(0, "cell_type", cell_type)
    if df.shape[1] > 7:
        df = df.iloc[:, :7]
    merged.append(df)

result = pd.concat(merged, ignore_index=True)
out_path = Path("DEG_celltype_merged.csv")
result.to_csv(out_path, index=False)

print(f"Merge {len(merged)} DEG files，write in {out_path}.")

deg_summary = result["cell_type"].value_counts().to_dict()
adata_deg = ad.AnnData(
    X=None,
    obs=result,
    var=pd.DataFrame(index=[]),
    uns={
        "cell_types": sorted(result["cell_type"].unique()),
        "counts_per_cell_type": deg_summary,
    },
)
adata_deg.write("DEG_celltype_merged.h5ad", compression="gzip")
print("同时写出 DEG_celltype_merged.h5ad")
