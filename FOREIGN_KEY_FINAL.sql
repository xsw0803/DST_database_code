-- Add foreign key, to connect each table.
ALTER TABLE sample_metadata
ADD FOREIGN KEY (region_ID) REFERENCES brain_region(region_ID);

ALTER TABLE cell_metadata
ADD CONSTRAINT fk_sample_ID
FOREIGN KEY (sample_ID) REFERENCES sample_metadata(sample_ID);

ALTER TABLE pathway
	ADD CONSTRAINT fk_pathway_ligand
		FOREIGN KEY (ligand) REFERENCES secreted_receptor_protein_dictionary(gene_name),
	ADD CONSTRAINT fk_pathway_receptor
		FOREIGN KEY (receptor) REFERENCES secreted_receptor_protein_dictionary(gene_name);
        
ALTER TABLE regional_gene
ADD CONSTRAINT fk_gene_name
FOREIGN KEY (gene_name) REFERENCES secreted_receptor_protein_dictionary(gene_name);

ALTER TABLE sample_metadata
ADD CONSTRAINT fk_region_ID
FOREIGN KEY (region_ID) REFERENCES count_source(region_ID);

ALTER TABLE regional_gene
ADD CONSTRAINT fk_region_ID2
FOREIGN KEY (region_ID) REFERENCES count_source(region_ID);

ALTER TABLE hvg_EC
ADD CONSTRAINT fk_gene_index_ec
FOREIGN KEY (gene_index) REFERENCES regional_gene(gene_index);

ALTER TABLE hvg_MTG
ADD CONSTRAINT fk_gene_index_mtg
FOREIGN KEY (gene_index) REFERENCES regional_gene(gene_index);

ALTER TABLE hvg_PFC
ADD CONSTRAINT fk_gene_index_pfc
FOREIGN KEY (gene_index) REFERENCES regional_gene(gene_index);

ALTER TABLE deg_EC
ADD CONSTRAINT fk_gene_index_deg_ec
FOREIGN KEY (gene_index) REFERENCES regional_gene(gene_index);

ALTER TABLE deg_MTG
ADD CONSTRAINT fk_gene_index_mtg_ec
FOREIGN KEY (gene_index) REFERENCES regional_gene(gene_index);

ALTER TABLE deg_PFC
ADD CONSTRAINT fk_gene_index_pfc_ec
FOREIGN KEY (gene_index) REFERENCES regional_gene(gene_index);