import scanpy as sc

h5_path = "/Users/xsw0803/Desktop/Learn/DST2/ICA/Database/snRNA/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5"
adata = sc.read_10x_h5(h5_path)
adata.var_names_make_unique()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.write("normalized_expression.h5ad", compression="gzip")


