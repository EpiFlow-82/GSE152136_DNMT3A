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

# ============================================================
# Copy hg38 CpG reference bigBed into the hub
# ============================================================

cpg_bigbed_source <- here(
  "tracks",
  "reference",
  "DNMT3A_window_hg38_CpGs.bb"
)

cpg_bigbed_hub <- here(
  "hub",
  "data",
  "DNMT3A_window_hg38_CpGs.bb"
)

if (!file.exists(cpg_bigbed_source)) {
  stop(
    "CpG reference bigBed not found: ",
    cpg_bigbed_source
  )
}

cpg_copy_ok <- file.copy(
  from = cpg_bigbed_source,
  to = cpg_bigbed_hub,
  overwrite = TRUE
)

if (!cpg_copy_ok) {
  stop("Failed to copy CpG reference bigBed into hub/data")
}

cat("CpG reference bigBed copied: 1\n")
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
# Copy WGBS group-summary BigWigs into the hub
# ============================================================

wgbs_summary_source <- c(
  AML = here(
    "tracks",
    "group_summaries",
    "WGBS_AML_mean_n16of18_DNMT3A.bw"
  ),
  PBMC = here(
    "tracks",
    "group_summaries",
    "WGBS_PBMC_mean_n2of2_DNMT3A.bw"
  )
)

stopifnot(
  all(file.exists(wgbs_summary_source))
)

wgbs_summary_hub_files <- file.path(
  hub_data_dir,
  basename(wgbs_summary_source)
)

summary_copy_ok <- file.copy(
  from = wgbs_summary_source,
  to = wgbs_summary_hub_files,
  overwrite = TRUE
)

cat(
  "\nWGBS summary BigWigs copied:",
  sum(summary_copy_ok),
  "of",
  length(summary_copy_ok),
  "\n"
)

# ============================================================
# Copy ATAC group-summary and derived-locus files into the hub
# ============================================================

atac_summary_source <- c(
  AML = here(
    "tracks",
    "ATAC_group_summary",
    "ATAC_AML_mean_RPGC_DNMT3A.bw"
  ),
  PBMC = here(
    "tracks",
    "ATAC_group_summary",
    "ATAC_PBMC_mean_RPGC_DNMT3A.bw"
  )
)

atac_loci_source <- here(
  "tracks",
  "ATAC_group_summary",
  "ATAC_derived_501bp_loci_DNMT3A.bb"
)

stopifnot(
  all(file.exists(atac_summary_source)),
  file.exists(atac_loci_source)
)

atac_summary_hub_files <- file.path(
  hub_data_dir,
  basename(atac_summary_source)
)

atac_loci_hub_file <- file.path(
  hub_data_dir,
  basename(atac_loci_source)
)

atac_summary_copy_ok <- file.copy(
  from = atac_summary_source,
  to = atac_summary_hub_files,
  overwrite = TRUE
)

atac_loci_copy_ok <- file.copy(
  from = atac_loci_source,
  to = atac_loci_hub_file,
  overwrite = TRUE
)

cat(
  "\nATAC group-summary BigWigs copied:",
  sum(atac_summary_copy_ok),
  "of",
  length(atac_summary_copy_ok),
  "\n"
)

cat(
  "ATAC derived-locus bigBed copied:",
  as.integer(atac_loci_copy_ok),
  "of 1\n"
)

cat("\nATAC hub-file QC:\n")

print(
  tibble(
    file = c(
      basename(atac_summary_hub_files),
      basename(atac_loci_hub_file)
    ),
    exists = c(
      file.exists(atac_summary_hub_files),
      file.exists(atac_loci_hub_file)
    ),
    size_bytes = c(
      file.info(atac_summary_hub_files)$size,
      file.info(atac_loci_hub_file)$size
    )
  ),
  n = Inf,
  width = Inf
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
# Composite assay tracks with AML/PBMC subgroups
# ============================================================

parent_definitions <- c(
  "track ATAC",
  "compositeTrack on",
  "shortLabel ATAC-seq",
  "longLabel ATAC-seq chromatin accessibility around DNMT3A (RPGC-normalized GEO processed data)",
  "type bigWig",
  "subGroup1 group Biological_group AML=AML PBMC=PBMC",
  "subGroup2 level Track_level Summary=Summary Individual=Individual",
  "dimensions dimensionX=group dimensionY=level",
  "sortOrder level=+ group=+",
  "visibility hide",
  "",
  
  "track CTCF",
  "compositeTrack on",
  "shortLabel CTCF",
  "longLabel CTCF CUT&Tag signal around DNMT3A",
  "type bigWig",
  "subGroup1 group Biological_group AML=AML PBMC=PBMC",
  "dimensions dimensionX=group",
  "sortOrder group=+",
  "visibility hide",
  "",
  
  "track H3K27ac",
  "compositeTrack on",
  "shortLabel H3K27ac",
  "longLabel H3K27ac CUT&Tag signal around DNMT3A",
  "type bigWig",
  "subGroup1 group Biological_group AML=AML PBMC=PBMC",
  "dimensions dimensionX=group",
  "sortOrder group=+",
  "visibility hide",
  "",
  
  "track H3K27me3",
  "compositeTrack on",
  "shortLabel H3K27me3",
  "longLabel H3K27me3 CUT&Tag signal around DNMT3A",
  "type bigWig",
  "subGroup1 group Biological_group AML=AML PBMC=PBMC",
  "dimensions dimensionX=group",
  "sortOrder group=+",
  "visibility hide",
  "",
  
  "track WGBS",
  "compositeTrack on",
  "shortLabel WGBS",
  "longLabel WGBS CG methylation around DNMT3A (10-bp bins; GEO processed data)",
  "type bigWig",
  "subGroup1 group Biological_group AML=AML PBMC=PBMC",
  "subGroup2 level Track_level Summary=Summary Individual=Individual",
  "dimensions dimensionX=group dimensionY=level",
  "sortOrder level=+ group=+",
  "visibility hide",
  ""
)

make_track_definition <- function(i) {
  
  x <- ucsc_tracks[i, ]
  
  track_color <- case_when(
    x$group == "AML" ~ "178,34,34",
    x$group == "PBMC" ~ "30,90,180"
  )
  
  display_settings <- if (x$track_type == "WGBS") {
    c(
      "autoScale off",
      "viewLimits 0:1"
    )
  } else {
    c(
      "autoScale on"
    )
  }
  
  c(
    paste("track", x$ucsc_track),
    paste("parent", x$parent_track, "off"),
    paste("shortLabel", x$short_label),
    paste("longLabel", x$long_label),
    if (x$track_type %in% c("WGBS", "ATAC-seq")) {
      paste0("subGroups group=", x$group, " level=Individual")
    } else {
      paste0("subGroups group=", x$group)
    },
    "type bigWig",
    paste("bigDataUrl", x$bigDataUrl),
    paste("color", track_color),
    "visibility dense",
    display_settings,
    "maxHeightPixels 100:32:8",
    ""
  )
}
# ============================================================
# ATAC group-summary track definitions
# ============================================================

atac_summary_definitions <- c(
  "track ATAC_AML_summary",
  "parent ATAC off",
  "shortLabel AML mean n=13",
  "longLabel AML mean ATAC-seq RPGC signal — 13 AML patients",
  "subGroups group=AML level=Summary",
  "type bigWig",
  "bigDataUrl data/ATAC_AML_mean_RPGC_DNMT3A.bw",
  "color 178,34,34",
  "visibility full",
  "autoScale off",
  "viewLimits 0:30",
  "maxHeightPixels 100:40:8",
  "",
  
  "track ATAC_PBMC_summary",
  "parent ATAC off",
  "shortLabel PBMC mean n=3",
  "longLabel PBMC mean ATAC-seq RPGC signal — 3 PBMC samples",
  "subGroups group=PBMC level=Summary",
  "type bigWig",
  "bigDataUrl data/ATAC_PBMC_mean_RPGC_DNMT3A.bw",
  "color 30,90,180",
  "visibility full",
  "autoScale off",
  "viewLimits 0:30",
  "maxHeightPixels 100:40:8",
  ""
)
# ============================================================
# WGBS group-summary track definitions
# ============================================================

wgbs_summary_definitions <- c(
  "track WGBS_AML_summary",
  "parent WGBS off",
  "shortLabel AML mean >=16/18",
  "longLabel AML mean WGBS CG methylation — >=16 of 18 samples per 10-bp bin",
  "subGroups group=AML level=Summary",
  "type bigWig",
  "bigDataUrl data/WGBS_AML_mean_n16of18_DNMT3A.bw",
  "color 178,34,34",
  "visibility full",
  "autoScale off",
  "viewLimits 0:1",
  "maxHeightPixels 100:40:8",
  "",
  
  "track WGBS_PBMC_summary",
  "parent WGBS off",
  "shortLabel PBMC mean 2/2",
  "longLabel PBMC mean WGBS CG methylation — 2 of 2 samples per 10-bp bin",
  "subGroups group=PBMC level=Summary",
  "type bigWig",
  "bigDataUrl data/WGBS_PBMC_mean_n2of2_DNMT3A.bw",
  "color 30,90,180",
  "visibility full",
  "autoScale off",
  "viewLimits 0:1",
  "maxHeightPixels 100:40:8",
  ""
)
child_definitions <- unlist(
  lapply(
    seq_len(nrow(ucsc_tracks)),
    make_track_definition
  )
)

# ============================================================
# Derived ATAC locus annotation track
# ============================================================

atac_loci_definition <- c(
  "track ATAC_derived_loci",
  "shortLabel ATAC loci",
  "longLabel Derived non-overlapping 501-bp ATAC loci with AML/PBMC descriptive summaries",
  "type bigBed 9 +",
  "bigDataUrl data/ATAC_derived_501bp_loci_DNMT3A.bb",
  "visibility pack",
  "color 80,80,80",
  "itemRgb on",
  ""
)

cpg_reference_definition <- c(
  "",
  "track CpG_reference_hg38",
  "shortLabel hg38 CpGs",
  "longLabel hg38 reference CpG dinucleotides around DNMT3A",
  "type bigBed 4",
  "bigDataUrl data/DNMT3A_window_hg38_CpGs.bb",
  "visibility dense",
  "color 60,60,60",
  ""
)
trackdb_txt <- c(
  parent_definitions,
  atac_summary_definitions,
  atac_loci_definition,
  wgbs_summary_definitions,
  cpg_reference_definition,
  child_definitions
)
writeLines(
  trackdb_txt,
  here("hub", "trackDb.txt")
)

cat("\nCreated trackDb.txt\n")
cat("Individual BigWig tracks:", nrow(ucsc_tracks), "\n")
cat("ATAC group-summary tracks:", 2, "\n")
cat("Derived ATAC locus tracks:", 1, "\n")
cat("WGBS group-summary tracks: 2\n")
cat("Composite assay tracks: 5\n")
cat("CpG reference tracks: 1\n")
cat("Total lines in trackDb.txt:", length(trackdb_txt), "\n")

cat("\nFirst 35 lines of trackDb.txt:\n\n")

cat(
  head(
    readLines(here("hub", "trackDb.txt")),
    35
  ),
  sep = "\n"
)