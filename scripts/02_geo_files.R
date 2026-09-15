# ============================================================
# GSE152136 DNMT3A multi-omics project
# Step 2: Inventory processed GEO supplementary files
# ============================================================

library(tidyverse)
library(here)
library(GEOquery)

# Load our validated metadata
master_samples <- read_tsv(
  here("metadata", "master_samples.tsv"),
  show_col_types = FALSE
)

sample_manifest <- read_tsv(
  here("metadata", "sample_manifest.tsv"),
  show_col_types = FALSE
)

# Verify
dim(sample_manifest)

sample_manifest |>
  count(assay, group)

# ============================================================
# Query supplementary files for each GSM sample
# WITHOUT downloading them
# ============================================================

get_gsm_files <- function(gsm_id) {
  
  message("Checking ", gsm_id, "...")
  
  x <- tryCatch(
    getGEOSuppFiles(
      gsm_id,
      makeDirectory = FALSE,
      fetch_files = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(x) || nrow(x) == 0) {
    return(
      tibble(
        GSM = gsm_id,
        fname = NA_character_,
        url = NA_character_
      )
    )
  }
  
  x |>
    as_tibble() |>
    mutate(GSM = gsm_id) |>
    select(GSM, fname, url)
}

gsm_files <- purrr::map_dfr(
  sample_manifest$GSM,
  get_gsm_files
)

gsm_files <- gsm_files |>
  left_join(
    sample_manifest |>
      select(GSM, GSE, sample_id, group, assay, mark, title),
    by = "GSM"
  )

dim(gsm_files)

gsm_files |>
  select(
    GSM,
    sample_id,
    group,
    assay,
    mark,
    fname
  ) |>
  head(50)

# ============================================================
# SAVE COMPLETE GEO FILE INVENTORY
# ============================================================

write_tsv(
  gsm_files,
  here("metadata", "geo_supplementary_files.tsv")
)

message("Complete GEO supplementary-file inventory saved.")

# ============================================================
# SELECT PROCESSED GENOMIC TRACK FILES
# ============================================================

processed_tracks <- gsm_files |>
  filter(
    case_when(
      assay == "WGBS" ~
        str_detect(fname, regex("\\.bw$", ignore_case = TRUE)),
      
      assay == "CUT&Tag" ~
        str_detect(fname, regex("\\.bigwig$", ignore_case = TRUE)),
      
      assay == "ATAC-seq" ~
        str_detect(fname, regex("\\.bigwig$", ignore_case = TRUE)),
      
      assay == "Hi-C" ~
        str_detect(fname, regex("\\.mcool$", ignore_case = TRUE)),
      
      TRUE ~ FALSE
    )
  )

write_tsv(
  processed_tracks,
  here("metadata", "processed_tracks.tsv")
)

message("Processed-track manifest saved.")

processed_tracks |>
  count(assay, group, mark)

# ============================================================
# SPLIT 1D TRACKS AND Hi-C
# ============================================================

tracks_1d <- processed_tracks |>
  filter(
    assay %in% c(
      "WGBS",
      "CUT&Tag",
      "ATAC-seq"
    )
  )

tracks_hic <- processed_tracks |>
  filter(
    assay == "Hi-C"
  )

write_tsv(
  tracks_1d,
  here("metadata", "processed_tracks_1d.tsv")
)

write_tsv(
  tracks_hic,
  here("metadata", "processed_tracks_hic.tsv")
)

message("1D and Hi-C track manifests saved.")











