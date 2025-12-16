-- Import data: Directly import csv files into database.
CREATE TABLE `1gene_dictionary_human` (
  `gene_name` text,
  `protein_name` text,
  `gene_type` text,
  `function_CC` text,
  `cognitive_function` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE `2lr_pairs_human` (
  `pathway_name` text,
  `ligand` text,
  `receptor` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE `3brain_region_ad` (
  `region_ID` text,
  `region_name` text,
  `AD_phase_involvement` text,
  `location` text,
  `function` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE `4cell_table_allregions` (
  `Index` int DEFAULT NULL,
  `Barcode` text,
  `Sample_ID` text,
  `Cell_type` text,
  `Cell_subtype` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE `5sample_table_allregions` (
  `Sample_ID` text,
  `Diagnosis` text,
  `Region_ID` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE `6regional_genes` (
  `gene_name` text,
  `region_ID` text
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE `7.1hvg_ec` (
  `gene_name` text,
  `mean` double DEFAULT NULL,
  `variance` double DEFAULT NULL,
  `variance_norm` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE `7.2hvg_mtg` (
  `gene_name` text,
  `means` double DEFAULT NULL,
  `variances` double DEFAULT NULL,
  `variances_norm` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE `7.3hvg_pfc` (
  `gene_name` text,
  `mean` double DEFAULT NULL,
  `variance` double DEFAULT NULL,
  `variance_norm` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE `8.1up_down_ec` (
  `cell_type` text,
  `gene_name` text,
  `score` double DEFAULT NULL,
  `log_fold_change` double DEFAULT NULL,
  `p_value` double DEFAULT NULL,
  `p_value_adj` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE `8.2up_down_mtg` (
  `cell_type` text,
  `gene_name` text,
  `score` double DEFAULT NULL,
  `log_fold_change` double DEFAULT NULL,
  `p_value` double DEFAULT NULL,
  `p_value_adj` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

CREATE TABLE `8.3up_down_pfc` (
  `cell_type` text,
  `gene_name` text,
  `score` double DEFAULT NULL,
  `log_fold_change` double DEFAULT NULL,
  `p_value` double DEFAULT NULL,
  `p_value_adj` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;