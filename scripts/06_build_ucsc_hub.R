library(tidyverse)
library(here)

bigwig_manifest <- read_tsv(
  here("metadata", "regional_1d_bigwig_tracks.tsv"),
  show_col_types = FALSE
)

hub_dir <- here("hub")

dir.create(
  hub_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Hub directory:\n")
cat(hub_dir, "\n\n")

cat("BigWig tracks available:", nrow(bigwig_manifest), "\n")

# ============================================================
# Create hub.txt
# ============================================================

hub_txt <- c(
  "hub GSE152136_DNMT3A",
  "shortLabel GSE152136 DNMT3A",
  "longLabel GSE152136 AML multi-omic tracks around DNMT3A",
  "genomesFile genomes.txt",
  "email Florian.Wolff@plus.ac.at"
)

writeLines(
  hub_txt,
  here("hub", "hub.txt")
)
# ============================================================
# Create genomes.txt
# ============================================================

genomes_txt <- c(
  "genome hg38",
  "trackDb trackDb.txt"
)

writeLines(
  genomes_txt,
  here("hub", "genomes.txt")
)

cat("\nCreated hub files:\n")
cat(here("hub", "hub.txt"), "\n")
cat(here("hub", "genomes.txt"), "\n")

cat("\nhub.txt:\n")
cat(readLines(here("hub", "hub.txt")), sep = "\n")

cat("\n\ngenomes.txt:\n")
cat(readLines(here("hub", "genomes.txt")), sep = "\n")

# ============================================================
# Inspect manifest before building trackDb.txt
# ============================================================

cat("\nManifest columns:\n")
print(names(bigwig_manifest))

cat("\nTrack counts by type and group:\n")

bigwig_manifest |>
  count(track_type, group) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nFirst 10 tracks:\n")

bigwig_manifest |>
  select(
    GSM,
    sample_id,
    group,
    track_type,
    title,
    bigwig_file
  ) |>
  slice_head(n = 10) |>
  as.data.frame() |>
  print(row.names = FALSE)
# ============================================================
# Copy final BigWigs into the self-contained hub
# ============================================================

hub_data_dir <- here("hub", "data")

dir.create(
  hub_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

copy_results <- bigwig_manifest |>
  mutate(
    hub_filename = basename(bigwig_file),
    hub_file = file.path(
      hub_data_dir,
      hub_filename
    )
  )

copy_ok <- file.copy(
  from = copy_results$bigwig_file,
  to = copy_results$hub_file,
  overwrite = TRUE
)

copy_results <- copy_results |>
  mutate(
    copy_ok = copy_ok
  )

cat("\nBigWig files copied:", sum(copy_results$copy_ok), "\n")
cat("Copy failures:", sum(!copy_results$copy_ok), "\n")

cat(
  "Files now in hub/data:",
  length(list.files(hub_data_dir, pattern = "\\.bw$")),
  "\n"
)

cat(
  "Total hub BigWig size:",
  round(
    sum(file.info(copy_results$hub_file)$size) / 1024^2,
    2
  ),
  "MB\n"
)

# ============================================================
# Prepare UCSC track information
# ============================================================

ucsc_tracks <- copy_results |>
  mutate(
    ucsc_track = paste0(
      "t_",
      gsub("[^A-Za-z0-9_]", "_", GSM),
      "_",
      gsub("[^A-Za-z0-9_]", "_", track_type)
    ),
    
    parent_track = case_when(
      track_type == "ATAC-seq" ~ "ATAC",
      track_type == "CTCF" ~ "CTCF",
      track_type == "H3K27ac" ~ "H3K27ac",
      track_type == "H3K27me3" ~ "H3K27me3",
      track_type == "WGBS" ~ "WGBS"
    ),
    
    short_label = paste(
      group,
      sample_id
    ),
    
    long_label = case_when(
      track_type == "WGBS" ~ paste0(
        group, " ", sample_id,
        " — WGBS CG methylation (10-bp bins; GEO processed data)"
      ),
      
      TRUE ~ paste0(
        group, " ", sample_id,
        " — ", track_type,
        " — GSE152136"
      )
    ),
    
    bigDataUrl = paste0(
      "data/",
      hub_filename
    )
  )
# ============================================================
# Parent tracks
# ============================================================

parent_definitions <- c(
  "track ATAC",
  "superTrack on",
  "shortLabel ATAC-seq",
  "longLabel ATAC-seq chromatin accessibility around DNMT3A",
  "",
  "track CTCF",
  "superTrack on",
  "shortLabel CTCF",
  "longLabel CTCF CUT&Tag signal around DNMT3A",
  "",
  "track H3K27ac",
  "superTrack on",
  "shortLabel H3K27ac",
  "longLabel H3K27ac CUT&Tag signal around DNMT3A",
  "",
  "track H3K27me3",
  "superTrack on",
  "shortLabel H3K27me3",
  "longLabel H3K27me3 CUT&Tag signal around DNMT3A",
  "",
  "track WGBS",
  "superTrack on",
  "shortLabel WGBS",
  "longLabel WGBS CG methylation around DNMT3A (10-bp bins)",
  ""
)

make_track_definition <- function(i) {
  
  x <- ucsc_tracks[i, ]
  
  c(
    paste("track", x$ucsc_track),
    paste("parent", x$parent_track),
    paste("shortLabel", x$short_label),
    paste("longLabel", x$long_label),
    "type bigWig",
    paste("bigDataUrl", x$bigDataUrl),
    "visibility dense",
    "autoScale on",
    "maxHeightPixels 100:32:8",
    ""
  )
}

child_definitions <- unlist(
  lapply(
    seq_len(nrow(ucsc_tracks)),
    make_track_definition
  )
)

trackdb_txt <- c(
  parent_definitions,
  child_definitions
)

writeLines(
  trackdb_txt,
  here("hub", "trackDb.txt")
)

cat("\nCreated trackDb.txt\n")
cat("Individual BigWig tracks:", nrow(ucsc_tracks), "\n")
cat("Parent assay tracks: 5\n")
cat("Total lines in trackDb.txt:", length(trackdb_txt), "\n")

cat("\nFirst 35 lines of trackDb.txt:\n\n")

cat(
  head(
    readLines(here("hub", "trackDb.txt")),
    35
  ),
  sep = "\n"
)










