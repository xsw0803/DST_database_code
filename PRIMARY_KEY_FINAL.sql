-- Add PRIMARY KEY：Rename tables for querying. Modify data type and length to add primary key.
ALTER TABLE `DST_database`.`1gene_dictionary_human` 
RENAME TO  `DST_database`.`secreted_receptor_protein_dictionary`,
MODIFY gene_name char(255), 
MODIFY gene_type char(255),
ADD PRIMARY KEY (gene_name);

CREATE TABLE pathway AS
SELECT DISTINCT ligand, receptor, pathway_name
FROM 2lr_pairs_human;

ALTER TABLE pathway
MODIFY  ligand char(255),
MODIFY receptor char(255),
MODIFY pathway_name char(255),
ADD PRIMARY KEY (ligand, receptor);

DROP TABLE 2lr_pairs_human;

ALTER TABLE `DST_database`.`3brain_region_ad` 
RENAME TO  `DST_database`.`brain_region`,
CHANGE `function` region_function text;

ALTER TABLE brain_region
MODIFY region_ID char(255),
MODIFY region_name char(255),
MODIFY AD_phase_involvement char(255),
MODIFY location char(255),
MODIFY region_function char(255),
ADD PRIMARY KEY (region_ID);

ALTER TABLE `DST_database`.`4cell_table_allregions` 
RENAME TO  `DST_database`.`cell_metadata`,
CHANGE COLUMN `Index` `sample_index` INT,
CHANGE Barcode barcode text,
CHANGE Sample_ID sample_ID text,
CHANGE Cell_type cell_type text,
CHANGE Cell_subtype cell_subtype text;

ALTER TABLE cell_metadata
MODIFY sample_index INT(255),
MODIFY barcode char(255),
MODIFY sample_ID char(255),
MODIFY cell_type char(255),
MODIFY cell_subtype char(255),
ADD PRIMARY KEY (sample_index);

ALTER TABLE `DST_database`.`5sample_table_allregions` 
RENAME TO  `DST_database`.`sample_metadata`,
CHANGE Sample_ID sample_ID text,
CHANGE Region_ID region_ID text,
CHANGE Diagnosis diagnosis text;

ALTER TABLE sample_metadata
MODIFY sample_ID char(255),
MODIFY region_ID char(255),
MODIFY diagnosis char(255),
ADD PRIMARY KEY (sample_ID);

ALTER TABLE `DST_database`.`6regional_genes` 
RENAME TO  `DST_database`.`regional_gene`,
MODIFY region_ID char(255),
MODIFY gene_name char(255),
ADD PRIMARY KEY (gene_index);

ALTER TABLE `DST_database`.`7.1hvg_ec` 
RENAME TO  `DST_database`.`hvg_EC`,
MODIFY gene_name char(255),
ADD PRIMARY KEY (gene_index);

ALTER TABLE `DST_database`.`7.2hvg_mtg` 
RENAME TO  `DST_database`.`hvg_MTG`,
CHANGE means mean DOUBLE,
CHANGE variances variance DOUBLE,
CHANGE variances_norm variance_norm DOUBLE;

ALTER TABLE hvg_MTG
MODIFY gene_name char(255),
ADD PRIMARY KEY (gene_index);

ALTER TABLE `DST_database`.`7.3hvg_pfc` 
RENAME TO  `DST_database`.`hvg_PFC`, 
MODIFY gene_name char(255),
ADD PRIMARY KEY (gene_index);

ALTER TABLE `DST_database`.`8.1up_down_ec` 
RENAME TO  `DST_database`.`deg_EC`,
MODIFY gene_name char(255),
MODIFY cell_type char(255),
ADD PRIMARY KEY (gene_name, cell_type);

ALTER TABLE `DST_database`.`8.2up_down_mtg` 
RENAME TO  `DST_database`.`deg_MTG`,
MODIFY gene_name char(255),
MODIFY cell_type char(255),
ADD PRIMARY KEY (gene_name, cell_type);

ALTER TABLE `DST_database`.`8.3up_down_pfc` 
RENAME TO  `DST_database`.`deg_PFC`,
MODIFY gene_name char(255),
MODIFY cell_type char(255),
ADD PRIMARY KEY (gene_name, cell_type);

ALTER TABLE count_source
ADD PRIMARY KEY (region_ID);

