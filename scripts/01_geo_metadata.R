# ============================================================
# GSE152136 DNMT3A multi-omics project
# Step 1: GEO metadata
# ============================================================

library(tidyverse)
library(data.table)
library(here)
library(GEOquery)

# GEO SuperSeries
geo_superseries <- "GSE152136"

# Relevant subseries
geo_series <- tibble(
  accession = c(
    "GSE152096",
    "GSE152099",
    "GSE152132",
    "GSE152134",
    "GSE152135"
  ),
  assay = c(
    "RNA-seq",
    "WGBS",
    "CUT&Tag",
    "ATAC-seq",
    "Hi-C"
  )
)

geo_series

# ============================================================
# Download GEO metadata for each subseries
# ============================================================

get_geo_metadata <- function(gse_id, assay_name) {
  
  gse <- getGEO(
    gse_id,
    GSEMatrix = TRUE,
    getGPL = FALSE
  )
  
  # Most GEO series return one ExpressionSet.
  # If multiple platforms exist, we take the first for now.
  eset <- gse[[1]]
  
  meta <- pData(eset) |>
    as.data.frame() |>
    tibble::rownames_to_column("GSM") |>
    tibble::as_tibble() |>
    mutate(
      GSE = gse_id,
      assay = assay_name
    )
  
  meta
}

metadata_list <- purrr::map2(
  geo_series$accession,
  geo_series$assay,
  get_geo_metadata
)

all_metadata <- bind_rows(metadata_list)

dim(all_metadata)

# ============================================================
# MASTER BIOLOGICAL SAMPLE CLASSIFICATION
#
# Biological identity follows the study sample classification.
# Assay availability will be determined independently from
# the actual records/files deposited in GEO.
# ============================================================

aml_samples <- c(
  "270", "546", "629", "773", "018", "472", "027", "168",
  "838", "424", "496", "798", "933", "504", "1071", "1324",
  "1916", "1360", "566", "260", "990", "936", "1021", "868",
  "1953"
)

pbmc_samples <- c(
  "673", "103", "576"
)

hspc_samples <- c(
  "HSPC-341",
  "HSPC-050",
  "HSPC-213",
  "HSPC-822"
)

master_samples <- tibble(
  sample_id = c(
    aml_samples,
    pbmc_samples,
    hspc_samples
  ),
  group = c(
    rep("AML", length(aml_samples)),
    rep("PBMC", length(pbmc_samples)),
    rep("CD34+ HSPC", length(hspc_samples))
  )
)

master_samples


# ============================================================
# CREATE CLEAN GEO SAMPLE MANIFEST
# ============================================================

sample_manifest <- all_metadata |>
  mutate(
    # Extract biological sample ID from GEO title
    sample_id = case_when(
      str_detect(title, "^HSPC-") ~
        str_extract(title, "^HSPC-[0-9]+"),
      
      TRUE ~
        str_extract(title, "^[0-9]+")
    ),
    
    # Identify CUT&Tag target
    mark = case_when(
      str_detect(title, "H3K27ac")  ~ "H3K27ac",
      str_detect(title, "H3K27me3") ~ "H3K27me3",
      str_detect(title, "CTCF")      ~ "CTCF",
      str_detect(title, "IgG")       ~ "IgG",
      TRUE                           ~ NA_character_
    )
  ) |>
  
  # Biological identity comes from our master sample table,
  # NOT from GEO source_name_ch1
  left_join(
    master_samples,
    by = "sample_id"
  ) |>
  
  select(
    GSM,
    GSE,
    sample_id,
    group,
    assay,
    mark,
    title,
    source_name_ch1
  )

# ============================================================
# SAVE VALIDATED METADATA
# ============================================================

write_tsv(
  master_samples,
  here("metadata", "master_samples.tsv")
)

write_tsv(
  sample_manifest,
  here("metadata", "sample_manifest.tsv")
)

message("Metadata files saved successfully.")

