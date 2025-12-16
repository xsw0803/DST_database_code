CREATE TABLE count_source
(region_ID char(255),
website char(255),
file_type char(255));

INSERT INTO count_source VALUES("EC", 
"https://huggingface.co/datasets/Glo-B/PC/resolve/main/preprocessed_normalized_counts.csv",
 "csv"),
("MTG",
"https://huggingface.co/datasets/Glo-B/PC/resolve/main/MTG_raw_gzip.h5ad",
"h5ad"),
("PFC",
"https://huggingface.co/datasets/Glo-B/PC/resolve/main/normalized_expression.h5ad",
"h5ad");


