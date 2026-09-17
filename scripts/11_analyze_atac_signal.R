# ============================================================
# 11_analyze_atac_signal.R
#
# Characterize regional ATAC-seq signal around DNMT3A
# before constructing AML/PBMC group summaries.
# ============================================================

library(tidyverse)
library(here)
library(cpp11bigwig)
library(GenomicRanges)
# ------------------------------------------------------------
# 1. Load regional 1D track manifest
# ------------------------------------------------------------

tracks_1d <- read_tsv(
  here("metadata", "regional_1d_bigwig_tracks.tsv"),
  show_col_types = FALSE
)

cat("Rows in regional 1D manifest:", nrow(tracks_1d), "\n")
cat("Columns:\n")
print(names(tracks_1d))

# Keep ATAC-seq tracks only
atac_tracks <- tracks_1d |>
  filter(assay == "ATAC-seq")

cat("\nATAC-seq tracks:", nrow(atac_tracks), "\n")

cat("\nATAC-seq tracks by biological group:\n")
print(
  atac_tracks |>
    count(group)
)

cat("\nATAC-seq samples:\n")
print(
  atac_tracks |>
    select(sample_id, group, everything())
)
cat("\nAssay values in manifest:\n")

print(
  tracks_1d |>
    count(assay, sort = TRUE)
)
# ------------------------------------------------------------
# 2. Read and characterize regional ATAC-seq BigWigs
# ------------------------------------------------------------

read_atac_bigwig <- function(file, sample_id, group) {
  
  bw <- cpp11bigwig::read_bigwig(file)
  
  tibble(
    sample_id = sample_id,
    group = group,
    chrom = bw$chrom,
    start = bw$start,
    end = bw$end,
    value = bw$value
  )
}

atac_signal <- purrr::pmap_dfr(
  atac_tracks |>
    select(bigwig_file, sample_id, group),
  \(bigwig_file, sample_id, group) {
    read_atac_bigwig(
      file = bigwig_file,
      sample_id = sample_id,
      group = group
    )
  }
)

cat("\nTotal stored ATAC intervals:",
    nrow(atac_signal), "\n")

cat("\nSignal value range:\n")
print(
  range(
    atac_signal$value,
    na.rm = TRUE
  )
)

cat("\nStored interval widths:\n")
print(
  atac_signal |>
    mutate(width = end - start) |>
    count(width, sort = TRUE) |>
    slice_head(n = 20)
)

cat("\nCoordinate modulo 10:\n")
print(
  atac_signal |>
    summarise(
      start_mod10_min = min(start %% 10),
      start_mod10_max = max(start %% 10),
      end_mod10_min = min(end %% 10),
      end_mod10_max = max(end %% 10)
    )
)

cat("\nPer-sample ATAC signal summary:\n")

atac_sample_summary <- atac_signal |>
  group_by(group, sample_id) |>
  summarise(
    n_intervals = n(),
    min_signal = min(value, na.rm = TRUE),
    median_signal = median(value, na.rm = TRUE),
    mean_signal = mean(value, na.rm = TRUE),
    p95_signal = quantile(value, 0.95, na.rm = TRUE),
    p99_signal = quantile(value, 0.99, na.rm = TRUE),
    max_signal = max(value, na.rm = TRUE),
    .groups = "drop"
  )

print(
  atac_sample_summary,
  n = Inf,
  width = Inf
)
# ------------------------------------------------------------
# 3. Inspect ATAC-seq coordinate grid and boundaries
# ------------------------------------------------------------

cat("\nPer-sample coordinate structure:\n")

atac_coordinate_summary <- atac_signal |>
  group_by(group, sample_id) |>
  summarise(
    first_start = min(start),
    last_end = max(end),
    n_intervals = n(),
    min_width = min(end - start),
    max_width = max(end - start),
    all_widths_multiple_10 = all((end - start) %% 10 == 0),
    start_mod10_values = paste(
      sort(unique(start %% 10)),
      collapse = ","
    ),
    end_mod10_values = paste(
      sort(unique(end %% 10)),
      collapse = ","
    ),
    .groups = "drop"
  )

print(
  atac_coordinate_summary,
  n = Inf,
  width = Inf
)

cat("\nIntervals whose start is not divisible by 10:\n")

atac_non10_starts <- atac_signal |>
  filter(start %% 10 != 0) |>
  mutate(
    width = end - start,
    start_mod10 = start %% 10,
    end_mod10 = end %% 10
  )

print(
  atac_non10_starts,
  n = 100,
  width = Inf
)

cat(
  "\nNumber of intervals with non-10-aligned starts:",
  nrow(atac_non10_starts),
  "\n"
)
# ------------------------------------------------------------
# 4. Reconstruct underlying 10-bp ATAC-seq signal grid
# ------------------------------------------------------------

atac_grid_start <- 24727880L
atac_grid_end   <- 25842590L

expand_atac_interval <- function(
    sample_id,
    group,
    chrom,
    start,
    end,
    value
) {
  
  # Restrict to the clean 10-bp-aligned analysis window
  clean_start <- max(start, atac_grid_start)
  clean_end   <- min(end, atac_grid_end)
  
  if (clean_start >= clean_end) {
    return(NULL)
  }
  
  bin_starts <- seq(
    from = clean_start,
    to = clean_end - 10L,
    by = 10L
  )
  
  tibble(
    sample_id = sample_id,
    group = group,
    chrom = chrom,
    start = bin_starts,
    end = bin_starts + 10L,
    value = value
  )
}

atac_10bp <- purrr::pmap_dfr(
  atac_signal,
  expand_atac_interval
)

cat("\nReconstructed ATAC 10-bp observations:",
    nrow(atac_10bp), "\n")

cat("\nCheck reconstructed bin widths:\n")
print(
  atac_10bp |>
    mutate(width = end - start) |>
    count(width)
)

cat("\nCheck coordinate alignment:\n")
print(
  atac_10bp |>
    summarise(
      start_mod10_min = min(start %% 10),
      start_mod10_max = max(start %% 10),
      end_mod10_min = min(end %% 10),
      end_mod10_max = max(end %% 10)
    )
)

cat("\nReconstructed bins per sample:\n")
print(
  atac_10bp |>
    count(group, sample_id, name = "n_bins"),
  n = Inf
)

cat("\nDuplicate sample/bin combinations:\n")

atac_duplicates <- atac_10bp |>
  count(sample_id, chrom, start, end) |>
  filter(n > 1)

cat(nrow(atac_duplicates), "\n")
# ------------------------------------------------------------
# 5. Calculate per-sample ATAC signal distributions
#    on the reconstructed common 10-bp grid
# ------------------------------------------------------------

atac_sample_10bp_summary <- atac_10bp |>
  group_by(group, sample_id) |>
  summarise(
    n_bins = n(),
    n_zero = sum(value == 0),
    pct_zero = 100 * mean(value == 0),
    min_signal = min(value),
    median_signal = median(value),
    mean_signal = mean(value),
    p75_signal = quantile(value, 0.75),
    p90_signal = quantile(value, 0.90),
    p95_signal = quantile(value, 0.95),
    p99_signal = quantile(value, 0.99),
    max_signal = max(value),
    .groups = "drop"
  )

cat("\nPer-sample ATAC signal summary on common 10-bp grid:\n")

print(
  atac_sample_10bp_summary,
  n = Inf,
  width = Inf
)

cat("\nRange of per-sample mean signal:\n")
print(
  range(atac_sample_10bp_summary$mean_signal)
)

cat("\nRange of per-sample 95th percentile:\n")
print(
  range(atac_sample_10bp_summary$p95_signal)
)

cat("\nRange of percentage zero-signal bins:\n")
print(
  range(atac_sample_10bp_summary$pct_zero)
)
write_tsv(
  atac_sample_10bp_summary,
  here(
    "metadata",
    "ATAC_DNMT3A_window_sample_signal_summary.tsv"
  )
)
# ------------------------------------------------------------
# 6. Pairwise correlation of ATAC-seq signal profiles
#    across the common 10-bp grid
# ------------------------------------------------------------

# Convert the long sample x bin table into a matrix:
# rows = genomic 10-bp bins
# columns = samples

atac_wide <- atac_10bp |>
  select(
    chrom,
    start,
    end,
    sample_id,
    value
  ) |>
  pivot_wider(
    names_from = sample_id,
    values_from = value
  ) |>
  arrange(chrom, start)

cat("\nDimensions of ATAC signal matrix:\n")
print(dim(atac_wide))

cat("\nMissing values in ATAC signal matrix:\n")
print(
  sum(is.na(atac_wide))
)

# Extract numeric signal matrix
atac_matrix <- atac_wide |>
  select(-chrom, -start, -end) |>
  as.matrix()

# Pearson correlation:
# similarity of quantitative signal profiles
atac_cor_pearson <- cor(
  atac_matrix,
  method = "pearson"
)

# Spearman correlation:
# similarity of relative/ranked signal profiles,
# less sensitive to very strong individual peaks
atac_cor_spearman <- cor(
  atac_matrix,
  method = "spearman"
)

cat("\nPearson correlation matrix:\n")
print(
  round(atac_cor_pearson, 3)
)

cat("\nSpearman correlation matrix:\n")
print(
  round(atac_cor_spearman, 3)
)

cat("\nRange of off-diagonal Pearson correlations:\n")
print(
  range(
    atac_cor_pearson[
      upper.tri(atac_cor_pearson)
    ]
  )
)

cat("\nRange of off-diagonal Spearman correlations:\n")
print(
  range(
    atac_cor_spearman[
      upper.tri(atac_cor_spearman)
    ]
  )
)
# ------------------------------------------------------------
# 7. Inspect deposited ATAC-seq supplementary files
# ------------------------------------------------------------

geo_files <- read_tsv(
  here("metadata", "geo_supplementary_files.tsv"),
  show_col_types = FALSE
)

cat("\nColumns in GEO supplementary-file inventory:\n")
print(names(geo_files))

cat("\nNumber of supplementary-file records:\n")
print(nrow(geo_files))

cat("\nATAC-related supplementary files:\n")

atac_geo_files <- geo_files |>
  filter(
    GSM %in% atac_tracks$GSM
  )

print(
  atac_geo_files,
  n = Inf,
  width = Inf
)
# ------------------------------------------------------------
# 8. Download and inspect one ATAC-seq peak file
# ------------------------------------------------------------

atac_peak_dir <- "C:/Genomics/GSE152136/ATAC_peaks"

dir.create(
  atac_peak_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Select PBMC 103 as pilot
peak_103 <- atac_geo_files |>
  filter(
    sample_id == "103",
    str_detect(fname, "MACS2\\.filt\\.peaks\\.txt\\.gz$")
  )

cat("\nPBMC 103 peak file:\n")
print(
  peak_103 |>
    select(GSM, sample_id, group, fname, url),
  width = Inf
)

peak_103_file <- file.path(
  atac_peak_dir,
  peak_103$fname
)

if (!file.exists(peak_103_file)) {
  
  download.file(
    url = peak_103$url,
    destfile = peak_103_file,
    mode = "wb",
    quiet = FALSE
  )
}

cat(
  "\nDownloaded peak file:\n",
  peak_103_file,
  "\n"
)

cat(
  "File size:",
  file.info(peak_103_file)$size,
  "bytes\n"
)

# Read first lines without assuming the column structure
peak_103_preview <- readLines(
  gzfile(peak_103_file),
  n = 10
)

cat("\nFirst 10 lines of PBMC 103 peak file:\n")
cat(
  peak_103_preview,
  sep = "\n"
)
# ------------------------------------------------------------
# 9. Download all 16 ATAC-seq MACS2 filtered peak files
# ------------------------------------------------------------

atac_peak_files <- atac_geo_files |>
  filter(
    str_detect(fname, "MACS2\\.filt\\.peaks\\.txt\\.gz$")
  ) |>
  arrange(group, sample_id)

cat("\nNumber of ATAC peak files expected:\n")
print(nrow(atac_peak_files))

cat("\nPeak files by biological group:\n")
print(
  atac_peak_files |>
    count(group)
)

download_atac_peak <- function(fname, url) {
  
  dest <- file.path(
    atac_peak_dir,
    fname
  )
  
  if (file.exists(dest)) {
    return("already_exists")
  }
  
  download.file(
    url = url,
    destfile = dest,
    mode = "wb",
    quiet = TRUE
  )
  
  if (file.exists(dest)) {
    return("success")
  }
  
  return("failed")
}

atac_peak_download_status <- atac_peak_files |>
  mutate(
    download_status = map2_chr(
      fname,
      url,
      download_atac_peak
    )
  )

cat("\nDownload status:\n")
print(
  atac_peak_download_status |>
    count(download_status)
)

cat("\nDownloaded files by group:\n")
print(
  atac_peak_download_status |>
    count(group, download_status)
)
# ------------------------------------------------------------
# 10. Parse and QC all ATAC-seq MACS2 filtered peak files
# ------------------------------------------------------------

read_atac_peaks <- function(fname, sample_id, group, GSM) {
  
  file_path <- file.path(
    atac_peak_dir,
    fname
  )
  
  read_tsv(
    file_path,
    col_names = c(
      "chrom",
      "start",
      "end",
      "peak_name",
      "score",
      "strand",
      "signal_value",
      "minus_log10_pvalue",
      "minus_log10_qvalue",
      "summit_offset"
    ),
    col_types = cols(
      chrom = col_character(),
      start = col_double(),
      end = col_double(),
      peak_name = col_character(),
      score = col_double(),
      strand = col_character(),
      signal_value = col_double(),
      minus_log10_pvalue = col_double(),
      minus_log10_qvalue = col_double(),
      summit_offset = col_double()
    ),
    show_col_types = FALSE
  ) |>
    mutate(
      GSM = GSM,
      sample_id = sample_id,
      group = group,
      peak_width = end - start,
      summit = start + summit_offset
    )
}

atac_peaks <- pmap_dfr(
  atac_peak_files |>
    select(fname, sample_id, group, GSM),
  read_atac_peaks
)

cat("\nTotal MACS2 filtered peaks across all 16 samples:\n")
print(nrow(atac_peaks))

cat("\nNumber of columns after parsing:\n")
print(ncol(atac_peaks))

cat("\nChromosome naming examples:\n")
print(
  sort(unique(atac_peaks$chrom))[1:10]
)

# Basic structural checks
cat("\nInvalid peak widths (end <= start):\n")
print(
  sum(atac_peaks$end <= atac_peaks$start)
)

cat("\nMissing values in essential peak columns:\n")
print(
  atac_peaks |>
    summarise(
      chrom_NA = sum(is.na(chrom)),
      start_NA = sum(is.na(start)),
      end_NA = sum(is.na(end)),
      signal_NA = sum(is.na(signal_value)),
      qvalue_NA = sum(is.na(minus_log10_qvalue)),
      summit_NA = sum(is.na(summit_offset))
    )
)

# Per-sample genome-wide peak QC
atac_peak_sample_summary <- atac_peaks |>
  group_by(group, sample_id, GSM) |>
  summarise(
    n_peaks = n(),
    median_width = median(peak_width),
    mean_width = mean(peak_width),
    p95_width = quantile(peak_width, 0.95),
    median_signal = median(signal_value),
    p95_signal = quantile(signal_value, 0.95),
    max_signal = max(signal_value),
    .groups = "drop"
  ) |>
  arrange(group, sample_id)

cat("\nGenome-wide peak summary by sample:\n")
print(
  atac_peak_sample_summary,
  n = Inf,
  width = Inf
)

# DNMT3A +/- 500 kb project window
dnmt3a_window_start <- 24727873L
dnmt3a_window_end   <- 25842590L

atac_peaks_dnmt3a <- atac_peaks |>
  filter(
    chrom == "chr2",
    start < dnmt3a_window_end,
    end > dnmt3a_window_start
  )

atac_dnmt3a_peak_counts <- atac_peaks_dnmt3a |>
  count(
    group,
    sample_id,
    GSM,
    name = "n_peaks_DNMT3A_window"
  ) |>
  right_join(
    atac_peak_files |>
      select(group, sample_id, GSM),
    by = c("group", "sample_id", "GSM")
  ) |>
  mutate(
    n_peaks_DNMT3A_window =
      replace_na(n_peaks_DNMT3A_window, 0L)
  ) |>
  arrange(group, sample_id)

cat("\nPeaks overlapping DNMT3A +/- 500 kb window:\n")
print(
  atac_dnmt3a_peak_counts,
  n = Inf,
  width = Inf
)

cat("\nTotal sample-specific peak calls in DNMT3A window:\n")
print(nrow(atac_peaks_dnmt3a))

# ------------------------------------------------------------
# 11. Identify recurrent ATAC peak regions in DNMT3A window
# ------------------------------------------------------------

# Convert DNMT3A-window peaks to GenomicRanges
atac_dnmt3a_gr <- GRanges(
  seqnames = atac_peaks_dnmt3a$chrom,
  ranges = IRanges(
    start = atac_peaks_dnmt3a$start + 1,
    end = atac_peaks_dnmt3a$end
  )
)

# Merge all overlapping peak intervals, regardless of sample
atac_consensus_gr <- reduce(atac_dnmt3a_gr)

cat("\nSample-specific peak calls:\n")
print(length(atac_dnmt3a_gr))

cat("\nMerged recurrent genomic regions:\n")
print(length(atac_consensus_gr))

cat("\nMerged-region width summary (bp):\n")
print(summary(width(atac_consensus_gr)))

# Determine which sample peaks overlap each merged region
peak_hits <- findOverlaps(
  atac_consensus_gr,
  atac_dnmt3a_gr
)

peak_support_long <- tibble(
  region_id = queryHits(peak_hits),
  peak_row = subjectHits(peak_hits)
) |>
  mutate(
    sample_id = atac_peaks_dnmt3a$sample_id[peak_row],
    group = atac_peaks_dnmt3a$group[peak_row]
  ) |>
  distinct(
    region_id,
    sample_id,
    group
  )

# Count independent samples supporting each region
peak_support <- peak_support_long |>
  count(
    region_id,
    group,
    name = "n_samples"
  ) |>
  pivot_wider(
    names_from = group,
    values_from = n_samples,
    values_fill = 0
  )

# Make sure both group columns exist
if (!"AML" %in% names(peak_support)) {
  peak_support$AML <- 0L
}

if (!"PBMC" %in% names(peak_support)) {
  peak_support$PBMC <- 0L
}

peak_support <- peak_support |>
  mutate(
    chrom = as.character(seqnames(atac_consensus_gr))[region_id],
    start = start(atac_consensus_gr)[region_id] - 1,
    end = end(atac_consensus_gr)[region_id],
    width = end - start,
    AML_fraction = AML / 13,
    PBMC_fraction = PBMC / 3
  ) |>
  select(
    region_id,
    chrom,
    start,
    end,
    width,
    AML,
    AML_fraction,
    PBMC,
    PBMC_fraction
  ) |>
  arrange(start)

cat("\nFirst recurrent regions and sample support:\n")
print(
  peak_support,
  n = 20,
  width = Inf
)

cat("\nAML support-count distribution:\n")
print(
  peak_support |>
    count(AML) |>
    arrange(AML)
)

cat("\nPBMC support-count distribution:\n")
print(
  peak_support |>
    count(PBMC) |>
    arrange(PBMC)
)

cat("\nRegions supported by at least 50% of each group:\n")

print(
  peak_support |>
    summarise(
      AML_at_least_7_of_13 = sum(AML >= 7),
      PBMC_at_least_2_of_3 = sum(PBMC >= 2)
    )
)

# ------------------------------------------------------------
# 12. Examine summit coherence within recurrent ATAC regions
# ------------------------------------------------------------

# Add each individual MACS2 peak summit to the recurrent region
summit_support_long <- tibble(
  region_id = queryHits(peak_hits),
  peak_row = subjectHits(peak_hits)
) |>
  mutate(
    sample_id = atac_peaks_dnmt3a$sample_id[peak_row],
    group = atac_peaks_dnmt3a$group[peak_row],
    peak_start = atac_peaks_dnmt3a$start[peak_row],
    peak_end = atac_peaks_dnmt3a$end[peak_row],
    summit = atac_peaks_dnmt3a$summit[peak_row]
  ) |>
  distinct(
    region_id,
    sample_id,
    group,
    .keep_all = TRUE
  )

# Check whether any sample contributes more than one peak
# to the same recurrent region
summit_duplicate_check <- tibble(
  region_id = queryHits(peak_hits),
  peak_row = subjectHits(peak_hits)
) |>
  mutate(
    sample_id = atac_peaks_dnmt3a$sample_id[peak_row]
  ) |>
  count(
    region_id,
    sample_id,
    name = "n_peaks"
  ) |>
  filter(n_peaks > 1)

cat("\nSample/region combinations with >1 peak:\n")
print(nrow(summit_duplicate_check))

if (nrow(summit_duplicate_check) > 0) {
  print(
    summit_duplicate_check,
    n = Inf
  )
}

# Summit spread within each recurrent region
summit_region_summary <- summit_support_long |>
  group_by(region_id) |>
  summarise(
    n_samples = n_distinct(sample_id),
    n_AML = n_distinct(sample_id[group == "AML"]),
    n_PBMC = n_distinct(sample_id[group == "PBMC"]),
    summit_min = min(summit),
    summit_max = max(summit),
    summit_span = summit_max - summit_min,
    summit_median = median(summit),
    .groups = "drop"
  ) |>
  mutate(
    chrom = as.character(seqnames(atac_consensus_gr))[region_id],
    region_start = start(atac_consensus_gr)[region_id] - 1,
    region_end = end(atac_consensus_gr)[region_id],
    region_width = region_end - region_start
  ) |>
  select(
    region_id,
    chrom,
    region_start,
    region_end,
    region_width,
    n_samples,
    n_AML,
    n_PBMC,
    summit_min,
    summit_max,
    summit_span,
    summit_median
  ) |>
  arrange(region_start)

cat("\nSummit-span summary across all 155 recurrent regions:\n")
print(
  summary(summit_region_summary$summit_span)
)

# Summit coherence is only informative when >=2 independent
# samples contribute peaks to the same region
summit_shared_summary <- summit_region_summary |>
  filter(n_samples >= 2)

cat("\nRegions supported by >=2 independent samples:\n")
print(nrow(summit_shared_summary))

cat("\nSummit-span summary for regions supported by >=2 samples:\n")
print(
  summary(summit_shared_summary$summit_span)
)

cat("\nMost widely separated summit regions:\n")
print(
  summit_shared_summary |>
    arrange(desc(summit_span)) |>
    slice_head(n = 15),
  n = 15,
  width = Inf
)

# ------------------------------------------------------------
# DNMT3A gene and exact intron 2
# ------------------------------------------------------------

dnmt3a_gene_start <- 25227873L
dnmt3a_gene_end   <- 25342590L

dnmt3a_intron2_start <- 25300243L
dnmt3a_intron2_end   <- 25313912L

summit_regions_gene <- summit_region_summary |>
  filter(
    region_start < dnmt3a_gene_end,
    region_end > dnmt3a_gene_start
  )

summit_regions_intron2 <- summit_region_summary |>
  filter(
    region_start < dnmt3a_intron2_end,
    region_end > dnmt3a_intron2_start
  )

cat("\nRecurrent regions overlapping DNMT3A gene:\n")
print(
  summit_regions_gene,
  n = Inf,
  width = Inf
)

cat("\nRecurrent regions overlapping exact DNMT3A intron 2:\n")
print(
  summit_regions_intron2,
  n = Inf,
  width = Inf
)
# ------------------------------------------------------------
# 12b. Corrected summit analysis: retain every MACS2 peak
# ------------------------------------------------------------

# Keep every peak that overlaps each union region.
# Do NOT collapse multiple peaks from the same patient.
summit_all_peaks <- tibble(
  region_id = queryHits(peak_hits),
  peak_row = subjectHits(peak_hits)
) |>
  mutate(
    sample_id = atac_peaks_dnmt3a$sample_id[peak_row],
    group = atac_peaks_dnmt3a$group[peak_row],
    peak_start = atac_peaks_dnmt3a$start[peak_row],
    peak_end = atac_peaks_dnmt3a$end[peak_row],
    peak_name = atac_peaks_dnmt3a$peak_name[peak_row],
    summit = atac_peaks_dnmt3a$summit[peak_row]
  )

# Correct region-level summary:
# n_peaks = all MACS2 peaks
# n_samples = distinct patients
summit_region_summary_corrected <- summit_all_peaks |>
  group_by(region_id) |>
  summarise(
    n_peaks = n(),
    n_samples = n_distinct(sample_id),
    n_AML = n_distinct(sample_id[group == "AML"]),
    n_PBMC = n_distinct(sample_id[group == "PBMC"]),
    n_unique_summits = n_distinct(summit),
    summit_min = min(summit),
    summit_max = max(summit),
    summit_span = summit_max - summit_min,
    summit_median = median(summit),
    summit_IQR = IQR(summit),
    .groups = "drop"
  ) |>
  mutate(
    chrom = as.character(seqnames(atac_consensus_gr))[region_id],
    region_start = start(atac_consensus_gr)[region_id] - 1,
    region_end = end(atac_consensus_gr)[region_id],
    region_width = region_end - region_start,
    summit_span_fraction = summit_span / region_width
  ) |>
  select(
    region_id,
    chrom,
    region_start,
    region_end,
    region_width,
    n_peaks,
    n_samples,
    n_AML,
    n_PBMC,
    n_unique_summits,
    summit_min,
    summit_max,
    summit_span,
    summit_IQR,
    summit_median,
    summit_span_fraction
  ) |>
  arrange(region_start)

cat("\nCorrected summit-span summary, all 155 union regions:\n")
print(
  summary(summit_region_summary_corrected$summit_span)
)

cat("\nCorrected summit-span summary, regions with >=2 patients:\n")
print(
  summit_region_summary_corrected |>
    filter(n_samples >= 2) |>
    pull(summit_span) |>
    summary()
)

cat("\nRegions containing more peaks than supporting patients:\n")
print(
  summit_region_summary_corrected |>
    filter(n_peaks > n_samples) |>
    arrange(desc(n_peaks - n_samples)),
  n = Inf,
  width = Inf
)

# Inspect all actual peaks/summits in especially broad or complex regions
regions_to_inspect <- c(7L, 89L, 99L, 137L)

cat("\nAll peaks and summits in regions 7, 89, 99 and 137:\n")
print(
  summit_all_peaks |>
    filter(region_id %in% regions_to_inspect) |>
    arrange(region_id, summit, sample_id) |>
    select(
      region_id,
      sample_id,
      group,
      peak_name,
      peak_start,
      peak_end,
      summit
    ),
  n = Inf,
  width = Inf
)

# Corrected summaries specifically within DNMT3A
summit_gene_corrected <- summit_region_summary_corrected |>
  filter(
    region_start < dnmt3a_gene_end,
    region_end > dnmt3a_gene_start
  )

summit_intron2_corrected <- summit_region_summary_corrected |>
  filter(
    region_start < dnmt3a_intron2_end,
    region_end > dnmt3a_intron2_start
  )

cat("\nCorrected DNMT3A gene regions:\n")
print(
  summit_gene_corrected,
  n = Inf,
  width = Inf
)

cat("\nCorrected exact intron 2 regions:\n")
print(
  summit_intron2_corrected,
  n = Inf,
  width = Inf
)

cat("\nAll individual MACS2 peaks overlapping exact intron 2:\n")
print(
  summit_all_peaks |>
    filter(region_id %in% summit_intron2_corrected$region_id) |>
    arrange(region_id, summit, sample_id) |>
    select(
      region_id,
      sample_id,
      group,
      peak_name,
      peak_start,
      peak_end,
      summit
    ),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 13. AML vs PBMC quantitative ATAC signal across full window
# ------------------------------------------------------------

# Calculate group summaries at every common 10-bp bin
atac_group_10bp <- atac_10bp |>
  group_by(chrom, start, end, group) |>
  summarise(
    n_samples = n(),
    mean_signal = mean(value),
    median_signal = median(value),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = group,
    values_from = c(n_samples, mean_signal, median_signal)
  ) |>
  mutate(
    mean_diff_AML_minus_PBMC =
      mean_signal_AML - mean_signal_PBMC,
    
    median_diff_AML_minus_PBMC =
      median_signal_AML - median_signal_PBMC
  ) |>
  arrange(chrom, start)

cat("\nNumber of 10-bp bins in AML vs PBMC comparison:\n")
print(nrow(atac_group_10bp))

cat("\nGroup sample counts represented per bin:\n")
print(
  atac_group_10bp |>
    count(n_samples_AML, n_samples_PBMC)
)

cat("\nRange of AML group summaries:\n")
print(
  summary(
    atac_group_10bp |>
      select(mean_signal_AML, median_signal_AML)
  )
)

cat("\nRange of PBMC group summaries:\n")
print(
  summary(
    atac_group_10bp |>
      select(mean_signal_PBMC, median_signal_PBMC)
  )
)

cat("\nAML minus PBMC difference summaries:\n")
print(
  summary(
    atac_group_10bp |>
      select(
        mean_diff_AML_minus_PBMC,
        median_diff_AML_minus_PBMC
      )
  )
)

# ------------------------------------------------------------
# Inspect strongest quantitative differences across full window
# ------------------------------------------------------------

cat("\n20 bins with largest positive AML - PBMC mean difference:\n")
print(
  atac_group_10bp |>
    arrange(desc(mean_diff_AML_minus_PBMC)) |>
    slice_head(n = 20),
  n = 20,
  width = Inf
)

cat("\n20 bins with largest negative AML - PBMC mean difference:\n")
print(
  atac_group_10bp |>
    arrange(mean_diff_AML_minus_PBMC) |>
    slice_head(n = 20),
  n = 20,
  width = Inf
)

cat("\n20 bins with largest positive AML - PBMC median difference:\n")
print(
  atac_group_10bp |>
    arrange(desc(median_diff_AML_minus_PBMC)) |>
    slice_head(n = 20),
  n = 20,
  width = Inf
)

cat("\n20 bins with largest negative AML - PBMC median difference:\n")
print(
  atac_group_10bp |>
    arrange(median_diff_AML_minus_PBMC) |>
    slice_head(n = 20),
  n = 20,
  width = Inf
)

# ------------------------------------------------------------
# Focused subsets for interpretation
# ------------------------------------------------------------

atac_group_gene <- atac_group_10bp |>
  filter(
    start < dnmt3a_gene_end,
    end > dnmt3a_gene_start
  )

atac_group_intron2 <- atac_group_10bp |>
  filter(
    start < dnmt3a_intron2_end,
    end > dnmt3a_intron2_start
  )

cat("\n10-bp bins overlapping DNMT3A gene:\n")
print(nrow(atac_group_gene))

cat("\n10-bp bins overlapping exact intron 2:\n")
print(nrow(atac_group_intron2))

# ------------------------------------------------------------
# Quantitative signal specifically around intron-2 peak regions
# ------------------------------------------------------------

region90_start <- 25301395L
region90_end   <- 25302202L

region91_start <- 25309488L
region91_end   <- 25309651L

atac_region90 <- atac_group_10bp |>
  filter(
    start < region90_end,
    end > region90_start
  )

atac_region91 <- atac_group_10bp |>
  filter(
    start < region91_end,
    end > region91_start
  )

summarise_atac_region <- function(x, region_name) {
  
  x |>
    summarise(
      region = region_name,
      
      n_bins = n(),
      
      AML_mean_of_bins =
        mean(mean_signal_AML),
      
      PBMC_mean_of_bins =
        mean(mean_signal_PBMC),
      
      AML_median_of_bins =
        median(median_signal_AML),
      
      PBMC_median_of_bins =
        median(median_signal_PBMC),
      
      mean_diff =
        mean(mean_diff_AML_minus_PBMC),
      
      median_diff =
        median(median_diff_AML_minus_PBMC)
    )
}

intron2_peak_signal_summary <- bind_rows(
  summarise_atac_region(
    atac_region90,
    "Intron2_region90"
  ),
  summarise_atac_region(
    atac_region91,
    "Intron2_region91"
  )
)

cat("\nQuantitative ATAC signal at intron-2 peak regions:\n")
print(
  intron2_peak_signal_summary,
  width = Inf
)

# ------------------------------------------------------------
# Save whole-window quantitative comparison
# ------------------------------------------------------------

readr::write_tsv(
  atac_group_10bp,
  here::here(
    "metadata",
    "ATAC_AML_vs_PBMC_10bp_signal_summary.tsv"
  )
)

cat(
  "\nSaved:",
  here::here(
    "metadata",
    "ATAC_AML_vs_PBMC_10bp_signal_summary.tsv"
  ),
  "\n"
)

# ------------------------------------------------------------
# 14. Integrate quantitative ATAC signal with peak-support loci
# ------------------------------------------------------------

# First inspect the Step 12b summit-summary object
cat("\nColumns available in corrected summit-region summary:\n")
print(names(summit_region_summary_corrected))

# ------------------------------------------------------------
# 14a. Quantitative RPGC signal for each patient and union region
# ------------------------------------------------------------

# Assign each 10-bp ATAC bin to overlapping union regions
atac_locus_sample_signal <- summit_region_summary_corrected |>
  select(
    region_id,
    chrom,
    region_start,
    region_end,
    region_width
  ) |>
  tidyr::crossing(
    atac_10bp |>
      select(
        sample_id,
        group,
        chrom_bin = chrom,
        start,
        end,
        value
      )
  ) |>
  filter(
    chrom == chrom_bin,
    start < region_end,
    end > region_start
  ) |>
  group_by(
    region_id,
    chrom,
    region_start,
    region_end,
    region_width,
    sample_id,
    group
  ) |>
  summarise(
    n_bins = n(),
    mean_RPGC = mean(value),
    median_RPGC = median(value),
    max_RPGC = max(value),
    .groups = "drop"
  )

cat("\nPatient × locus quantitative observations:\n")
print(nrow(atac_locus_sample_signal))

cat("\nNumber of patients represented per locus:\n")
print(
  atac_locus_sample_signal |>
    count(region_id, name = "n_patients") |>
    count(n_patients)
)

# ------------------------------------------------------------
# 14b. AML and PBMC quantitative summaries per union region
# ------------------------------------------------------------

atac_locus_group_signal <- atac_locus_sample_signal |>
  group_by(
    region_id,
    chrom,
    region_start,
    region_end,
    region_width,
    group
  ) |>
  summarise(
    n_signal_samples = n_distinct(sample_id),
    
    mean_patient_RPGC =
      mean(mean_RPGC),
    
    median_patient_RPGC =
      median(mean_RPGC),
    
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = group,
    values_from = c(
      n_signal_samples,
      mean_patient_RPGC,
      median_patient_RPGC
    )
  ) |>
  mutate(
    mean_RPGC_diff_AML_minus_PBMC =
      mean_patient_RPGC_AML -
      mean_patient_RPGC_PBMC,
    
    median_RPGC_diff_AML_minus_PBMC =
      median_patient_RPGC_AML -
      median_patient_RPGC_PBMC
  )

cat("\nNumber of union regions with quantitative summaries:\n")
print(nrow(atac_locus_group_signal))

cat("\nSignal sample counts represented per locus:\n")
print(
  atac_locus_group_signal |>
    count(
      n_signal_samples_AML,
      n_signal_samples_PBMC
    )
)

# ------------------------------------------------------------
# 14c. Join quantitative signal with MACS2 patient support
# ------------------------------------------------------------

atac_locus_integrated <- summit_region_summary_corrected |>
  left_join(
    atac_locus_group_signal,
    by = c(
      "region_id",
      "chrom",
      "region_start",
      "region_end",
      "region_width"
    )
  ) |>
  mutate(
    AML_peak_fraction = n_AML / 13,
    PBMC_peak_fraction = n_PBMC / 3,
    
    peak_fraction_diff_AML_minus_PBMC =
      AML_peak_fraction - PBMC_peak_fraction
  ) |>
  arrange(region_id)

cat("\nIntegrated ATAC locus table:\n")
print(
  atac_locus_integrated |>
    select(
      region_id,
      chrom,
      region_start,
      region_end,
      region_width,
      n_peaks,
      n_samples,
      n_AML,
      n_PBMC,
      AML_peak_fraction,
      PBMC_peak_fraction,
      mean_patient_RPGC_AML,
      mean_patient_RPGC_PBMC,
      mean_RPGC_diff_AML_minus_PBMC,
      median_patient_RPGC_AML,
      median_patient_RPGC_PBMC,
      median_RPGC_diff_AML_minus_PBMC,
      n_unique_summits,
      summit_span
    ),
  n = 20,
  width = Inf
)

# ------------------------------------------------------------
# 14d. Inspect strongest quantitative differences
# ------------------------------------------------------------

cat("\n15 loci with largest positive AML - PBMC mean RPGC difference:\n")
print(
  atac_locus_integrated |>
    arrange(desc(mean_RPGC_diff_AML_minus_PBMC)) |>
    select(
      region_id,
      chrom,
      region_start,
      region_end,
      n_AML,
      n_PBMC,
      AML_peak_fraction,
      PBMC_peak_fraction,
      mean_patient_RPGC_AML,
      mean_patient_RPGC_PBMC,
      mean_RPGC_diff_AML_minus_PBMC,
      median_RPGC_diff_AML_minus_PBMC,
      n_peaks,
      n_samples,
      n_unique_summits,
      summit_span
    ) |>
    slice_head(n = 15),
  n = 15,
  width = Inf
)

cat("\n15 loci with largest negative AML - PBMC mean RPGC difference:\n")
print(
  atac_locus_integrated |>
    arrange(mean_RPGC_diff_AML_minus_PBMC) |>
    select(
      region_id,
      chrom,
      region_start,
      region_end,
      n_AML,
      n_PBMC,
      AML_peak_fraction,
      PBMC_peak_fraction,
      mean_patient_RPGC_AML,
      mean_patient_RPGC_PBMC,
      mean_RPGC_diff_AML_minus_PBMC,
      median_RPGC_diff_AML_minus_PBMC,
      n_peaks,
      n_samples,
      n_unique_summits,
      summit_span
    ) |>
    slice_head(n = 15),
  n = 15,
  width = Inf
)

# ------------------------------------------------------------
# 14e. Inspect DNMT3A gene and exact intron 2
# ------------------------------------------------------------

cat("\nIntegrated loci overlapping DNMT3A gene:\n")
print(
  atac_locus_integrated |>
    filter(
      region_start < dnmt3a_gene_end,
      region_end > dnmt3a_gene_start
    ) |>
    select(
      region_id,
      region_start,
      region_end,
      n_AML,
      n_PBMC,
      AML_peak_fraction,
      PBMC_peak_fraction,
      mean_patient_RPGC_AML,
      mean_patient_RPGC_PBMC,
      mean_RPGC_diff_AML_minus_PBMC,
      median_RPGC_diff_AML_minus_PBMC,
      n_peaks,
      n_samples,
      n_unique_summits,
      summit_span
    ),
  n = Inf,
  width = Inf
)

cat("\nIntegrated loci overlapping exact intron 2:\n")
print(
  atac_locus_integrated |>
    filter(
      region_start < dnmt3a_intron2_end,
      region_end > dnmt3a_intron2_start
    ) |>
    select(
      region_id,
      region_start,
      region_end,
      n_AML,
      n_PBMC,
      AML_peak_fraction,
      PBMC_peak_fraction,
      mean_patient_RPGC_AML,
      mean_patient_RPGC_PBMC,
      mean_RPGC_diff_AML_minus_PBMC,
      median_patient_RPGC_AML,
      median_patient_RPGC_PBMC,
      median_RPGC_diff_AML_minus_PBMC,
      n_peaks,
      n_samples,
      n_unique_summits,
      summit_span
    ),
  n = Inf,
  width = Inf
)
# ------------------------------------------------------------
# 15a. Inspect peak-level object before summit-aware refinement
# ------------------------------------------------------------

cat("\nObjects containing 'peak' in the current R session:\n")
print(ls(pattern = "[Pp][Ee][Aa][Kk]"))

cat("\nColumns in atac_peaks:\n")
print(names(atac_peaks))

cat("\nFirst 6 rows of atac_peaks:\n")
print(head(atac_peaks))

# ------------------------------------------------------------
# 15b. Compare summit-centered candidate locus widths
# ------------------------------------------------------------

# Keep only peaks overlapping the full DNMT3A +/-500 kb window
atac_peaks_window <- atac_peaks |>
  filter(
    chrom == "chr2",
    start < 25842590,
    end > 24727873
  )

cat("\nATAC peaks in full DNMT3A +/-500 kb window:\n")
print(nrow(atac_peaks_window))

cat("\nSamples represented:\n")
print(
  atac_peaks_window |>
    count(group, sample_id) |>
    count(group, name = "n_samples")
)

# Function:
# create fixed-width summit-centered intervals and reduce overlaps.
#
# Odd widths are used so the summit is the exact central base:
# 201 bp = summit +/-100 bp
# 301 bp = summit +/-150 bp
# 501 bp = summit +/-250 bp

build_summit_loci <- function(peaks, width_bp) {
  
  half_width <- (width_bp - 1L) / 2L
  
  summit_gr <- GenomicRanges::GRanges(
    seqnames = peaks$chrom,
    ranges = IRanges::IRanges(
      start = peaks$summit - half_width,
      end   = peaks$summit + half_width
    )
  )
  
  loci_gr <- GenomicRanges::reduce(
    summit_gr,
    ignore.strand = TRUE
  )
  
  tibble(
    chrom = as.character(GenomicRanges::seqnames(loci_gr)),
    start = GenomicRanges::start(loci_gr),
    end = GenomicRanges::end(loci_gr),
    width = GenomicRanges::width(loci_gr),
    tested_width = width_bp
  )
}

summit_loci_201 <- build_summit_loci(
  atac_peaks_window,
  201L
)

summit_loci_301 <- build_summit_loci(
  atac_peaks_window,
  301L
)

summit_loci_501 <- build_summit_loci(
  atac_peaks_window,
  501L
)

summit_width_comparison <- bind_rows(
  summit_loci_201,
  summit_loci_301,
  summit_loci_501
)

cat("\nNumber of resulting loci by tested summit width:\n")
print(
  summit_width_comparison |>
    count(tested_width, name = "n_loci")
)

cat("\nResulting locus-width distributions:\n")
print(
  summit_width_comparison |>
    group_by(tested_width) |>
    summarise(
      n_loci = n(),
      min_width = min(width),
      q25_width = quantile(width, 0.25),
      median_width = median(width),
      mean_width = mean(width),
      q75_width = quantile(width, 0.75),
      max_width = max(width),
      .groups = "drop"
    )
)

cat("\nNumber of loci wider than the nominal tested width:\n")
print(
  summit_width_comparison |>
    group_by(tested_width) |>
    summarise(
      n_loci = n(),
      n_wider_than_nominal =
        sum(width > tested_width),
      pct_wider_than_nominal =
        100 * mean(width > tested_width),
      .groups = "drop"
    )
)

# ------------------------------------------------------------
# 15c. Inspect summit-centered loci inside complex union regions
# ------------------------------------------------------------

complex_region_ids <- c(7L, 29L, 99L, 137L)

complex_regions <- summit_region_summary_corrected |>
  filter(region_id %in% complex_region_ids) |>
  select(
    region_id,
    chrom,
    region_start,
    region_end,
    region_width,
    n_peaks,
    n_samples,
    n_unique_summits,
    summit_min,
    summit_max,
    summit_span
  )

cat("\nOriginal complex union regions:\n")
print(complex_regions, n = Inf, width = Inf)

# Give candidate loci an ID within each tested width
summit_width_comparison_id <- summit_width_comparison |>
  group_by(tested_width) |>
  arrange(chrom, start, end, .by_group = TRUE) |>
  mutate(candidate_locus_id = row_number()) |>
  ungroup()


# Rename candidate-locus chromosome before crossing
summit_width_comparison_id2 <- summit_width_comparison_id |>
  dplyr::rename(candidate_chrom = chrom)

# Find candidate summit-centered loci overlapping each
# selected original union region
complex_candidate_loci <- tidyr::crossing(
  complex_regions,
  summit_width_comparison_id2
) |>
  filter(
    chrom == candidate_chrom,
    start <= region_end,
    end >= region_start
  ) |>
  select(
    region_id,
    chrom,
    original_start = region_start,
    original_end = region_end,
    original_width = region_width,
    original_n_peaks = n_peaks,
    original_n_samples = n_samples,
    original_n_unique_summits = n_unique_summits,
    original_summit_span = summit_span,
    tested_width,
    candidate_locus_id,
    candidate_start = start,
    candidate_end = end,
    candidate_width = width
  ) |>
  arrange(
    region_id,
    tested_width,
    candidate_start
  )

cat("\nCandidate loci overlapping complex regions:\n")
print(
  complex_candidate_loci,
  n = Inf,
  width = Inf
)

cat("\nNumber of candidate loci per complex region and tested width:\n")
print(
  complex_candidate_loci |>
    count(
      region_id,
      tested_width,
      name = "n_candidate_loci"
    ) |>
    arrange(
      region_id,
      tested_width
    ),
  n = Inf
)
# ------------------------------------------------------------
# 15d. Inspect empirical summit spacing in complex regions
# ------------------------------------------------------------

complex_peak_summits <- atac_peaks_window |>
  filter(
    chrom == "chr2"
  ) |>
  tidyr::crossing(
    complex_regions |>
      select(
        region_id,
        region_start,
        region_end
      )
  ) |>
  filter(
    start < region_end,
    end > region_start
  ) |>
  select(
    region_id,
    sample_id,
    group,
    peak_name,
    start,
    end,
    summit
  ) |>
  arrange(
    region_id,
    summit,
    sample_id
  )

# Calculate distance from each summit to the preceding
# distinct summit within the same original union region
complex_distinct_summits <- complex_peak_summits |>
  distinct(
    region_id,
    summit
  ) |>
  arrange(
    region_id,
    summit
  ) |>
  group_by(region_id) |>
  mutate(
    previous_summit = lag(summit),
    gap_to_previous_summit =
      summit - previous_summit
  ) |>
  ungroup()

cat("\nAll peak summits in the four complex regions:\n")
print(
  complex_peak_summits,
  n = Inf,
  width = Inf
)

cat("\nDistinct summit coordinates and consecutive gaps:\n")
print(
  complex_distinct_summits,
  n = Inf,
  width = Inf
)

cat("\nLargest consecutive summit gaps in each complex region:\n")
print(
  complex_distinct_summits |>
    filter(!is.na(gap_to_previous_summit)) |>
    group_by(region_id) |>
    arrange(
      desc(gap_to_previous_summit),
      .by_group = TRUE
    ) |>
    slice_head(n = 5) |>
    ungroup(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 16a. Global summit-gap distribution across the full window
# ------------------------------------------------------------

all_distinct_summits <- atac_peaks_window |>
  select(
    chrom,
    summit
  ) |>
  distinct() |>
  arrange(
    chrom,
    summit
  ) |>
  group_by(chrom) |>
  mutate(
    previous_summit = lag(summit),
    gap_to_previous_summit =
      summit - previous_summit
  ) |>
  ungroup()

all_summit_gaps <- all_distinct_summits |>
  filter(
    !is.na(gap_to_previous_summit)
  )

cat("\nNumber of peaks in regional dataset:\n")
print(nrow(atac_peaks_window))

cat("\nNumber of distinct summit coordinates:\n")
print(nrow(all_distinct_summits))

cat("\nGlobal consecutive-summit gap distribution (bp):\n")
print(
  summary(
    all_summit_gaps$gap_to_previous_summit
  )
)

cat("\nSelected summit-gap quantiles (bp):\n")
print(
  quantile(
    all_summit_gaps$gap_to_previous_summit,
    probs = c(
      0,
      0.10,
      0.25,
      0.50,
      0.75,
      0.90,
      0.95,
      0.99,
      1
    )
  )
)

cat("\nCounts of consecutive summit gaps by distance class:\n")

summit_gap_classes <- all_summit_gaps |>
  mutate(
    gap_class = cut(
      gap_to_previous_summit,
      breaks = c(
        -Inf,
        25,
        50,
        100,
        150,
        200,
        250,
        300,
        500,
        1000,
        Inf
      ),
      labels = c(
        "<=25",
        "26-50",
        "51-100",
        "101-150",
        "151-200",
        "201-250",
        "251-300",
        "301-500",
        "501-1000",
        ">1000"
      )
    )
  ) |>
  count(
    gap_class,
    .drop = FALSE,
    name = "n_gaps"
  ) |>
  mutate(
    pct_gaps =
      100 * n_gaps / sum(n_gaps)
  )

print(
  summit_gap_classes,
  n = Inf
)
# ------------------------------------------------------------
# 16b. ArchR-style iterative overlap removal
# ------------------------------------------------------------

# Rank MACS2 peak significance within each patient.
#
# minus_log10_qvalue is used because it is the deposited
# MACS2 q-value significance measure.
#
# percent_rank gives larger values to more significant peaks.

atac_peak_candidates <- atac_peaks_window |>
  group_by(sample_id) |>
  mutate(
    significance_rank =
      dplyr::percent_rank(minus_log10_qvalue)
  ) |>
  ungroup() |>
  arrange(
    desc(significance_rank),
    desc(minus_log10_qvalue),
    chrom,
    summit,
    sample_id
  ) |>
  mutate(
    candidate_id = row_number()
  )

cat("\nCandidate peak significance QC:\n")

print(
  atac_peak_candidates |>
    group_by(sample_id, group) |>
    summarise(
      n_peaks = n(),
      min_rank = min(significance_rank),
      median_rank = median(significance_rank),
      max_rank = max(significance_rank),
      .groups = "drop"
    ),
  n = Inf
)

# Create fixed 501-bp summit-centered candidate intervals.
#
# NOTE:
# summit coordinates in the GEO narrowPeak files are currently
# represented using the BED-derived absolute summit convention.
# We keep that convention consistently for this diagnostic step.
# Exact BED/GRanges coordinate conversion will be finalized before
# writing the UCSC track.

atac_peak_candidates <- atac_peak_candidates |>
  mutate(
    candidate_start = summit - 250L,
    candidate_end = summit + 250L
  )

# Iterative overlap removal:
# candidates are already ordered from highest to lowest
# normalized significance.
#
# Keep the best remaining candidate and discard all later
# candidates that overlap it.

candidate_gr <- GenomicRanges::GRanges(
  seqnames = atac_peak_candidates$chrom,
  ranges = IRanges::IRanges(
    start = atac_peak_candidates$candidate_start,
    end = atac_peak_candidates$candidate_end
  )
)

keep_candidate <- rep(FALSE, length(candidate_gr))
available <- rep(TRUE, length(candidate_gr))

for (i in seq_along(candidate_gr)) {
  
  if (!available[i]) {
    next
  }
  
  keep_candidate[i] <- TRUE
  
  overlapping <- GenomicRanges::findOverlaps(
    candidate_gr[i],
    candidate_gr,
    ignore.strand = TRUE
  )
  
  overlapping_indices <-
    S4Vectors::subjectHits(overlapping)
  
  available[overlapping_indices] <- FALSE
}

atac_archr_style_loci <- atac_peak_candidates |>
  filter(keep_candidate) |>
  arrange(
    chrom,
    candidate_start,
    candidate_end
  ) |>
  mutate(
    locus_id = row_number()
  )

cat("\nArchR-style iterative-overlap locus count:\n")
print(nrow(atac_archr_style_loci))

cat("\nWidth QC:\n")
print(
  atac_archr_style_loci |>
    summarise(
      n_loci = n(),
      min_width =
        min(candidate_end - candidate_start + 1L),
      median_width =
        median(candidate_end - candidate_start + 1L),
      max_width =
        max(candidate_end - candidate_start + 1L)
    )
)

# Confirm that retained loci do not overlap one another.

retained_gr <- GenomicRanges::GRanges(
  seqnames = atac_archr_style_loci$chrom,
  ranges = IRanges::IRanges(
    start = atac_archr_style_loci$candidate_start,
    end = atac_archr_style_loci$candidate_end
  )
)

retained_overlap_hits_all <- GenomicRanges::findOverlaps(
  retained_gr,
  retained_gr,
  ignore.strand = TRUE
)

retained_overlap_hits <- retained_overlap_hits_all[
  S4Vectors::queryHits(retained_overlap_hits_all) <
    S4Vectors::subjectHits(retained_overlap_hits_all)
]
cat("\nNumber of overlaps among retained loci:\n")
print(length(retained_overlap_hits))

# Inspect how many retained loci occur in the four
# previously problematic original union regions.

archr_complex_loci <- tidyr::crossing(
  complex_regions |>
    dplyr::rename(original_chrom = chrom),
  atac_archr_style_loci |>
    dplyr::rename(candidate_chrom = chrom)
) |>
  filter(
    original_chrom == candidate_chrom,
    candidate_start <= region_end,
    candidate_end >= region_start
  ) |>
  select(
    region_id,
    original_start = region_start,
    original_end = region_end,
    locus_id,
    candidate_start,
    candidate_end,
    summit,
    sample_id,
    group,
    minus_log10_qvalue,
    significance_rank
  ) |>
  arrange(
    region_id,
    candidate_start
  )

cat("\nRetained loci overlapping the four complex regions:\n")
print(
  archr_complex_loci,
  n = Inf,
  width = Inf
)

cat("\nNumber of retained loci per complex region:\n")
print(
  archr_complex_loci |>
    count(
      region_id,
      name = "n_retained_loci"
    ),
  n = Inf
)

# ------------------------------------------------------------
# 16c. Inspect MACS2 q-value ties before final peak ranking
# ------------------------------------------------------------

atac_qvalue_ties <- atac_peaks_window |>
  group_by(
    sample_id,
    group
  ) |>
  summarise(
    n_peaks = n(),
    n_unique_qvalues =
      n_distinct(minus_log10_qvalue),
    n_duplicated_qvalue_peaks =
      n_peaks - n_unique_qvalues,
    pct_duplicated_qvalue_peaks =
      100 * n_duplicated_qvalue_peaks / n_peaks,
    min_qvalue_score =
      min(minus_log10_qvalue),
    median_qvalue_score =
      median(minus_log10_qvalue),
    max_qvalue_score =
      max(minus_log10_qvalue),
    .groups = "drop"
  ) |>
  arrange(
    desc(pct_duplicated_qvalue_peaks)
  )

cat("\nMACS2 q-value tie summary by patient:\n")

print(
  atac_qvalue_ties,
  n = Inf,
  width = Inf
)

cat("\nMost frequent q-value scores within each patient:\n")

atac_qvalue_tie_detail <- atac_peaks_window |>
  count(
    sample_id,
    group,
    minus_log10_qvalue,
    name = "n_peaks"
  ) |>
  group_by(
    sample_id,
    group
  ) |>
  arrange(
    desc(n_peaks),
    desc(minus_log10_qvalue),
    .by_group = TRUE
  ) |>
  slice_head(n = 5) |>
  ungroup()

print(
  atac_qvalue_tie_detail,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 16d. Final deterministic fixed-width ATAC locus construction
# ------------------------------------------------------------

# Rank peaks within each patient.
#
# Primary criterion:
#   MACS2 -log10(q-value)
#
# Tie-breaker:
#   MACS2 signal_value
#
# Final deterministic tie-breaker:
#   genomic summit position
#
# The normalized rank is used ONLY to choose representative
# coordinates during overlap removal. It is not a biological
# locus score.

atac_peak_candidates_final <- atac_peaks_window |>
  group_by(sample_id) |>
  arrange(
    desc(minus_log10_qvalue),
    desc(signal_value),
    summit,
    .by_group = TRUE
  ) |>
  mutate(
    patient_peak_rank = row_number(),
    n_patient_peaks = n(),
    normalized_rank =
      if_else(
        n_patient_peaks > 1L,
        1 - (patient_peak_rank - 1) /
          (n_patient_peaks - 1),
        1
      )
  ) |>
  ungroup()

cat("\nFinal within-patient ranking QC:\n")

print(
  atac_peak_candidates_final |>
    group_by(sample_id, group) |>
    summarise(
      n_peaks = n(),
      best_rank = max(normalized_rank),
      median_rank = median(normalized_rank),
      worst_rank = min(normalized_rank),
      .groups = "drop"
    ),
  n = Inf
)

# Order all candidates for deterministic iterative selection.
#
# Across patients:
#   1. normalized within-patient rank
#   2. raw MACS2 q-value significance
#   3. MACS2 signal value
#   4. genomic summit
#   5. sample ID

atac_peak_candidates_final <- atac_peak_candidates_final |>
  arrange(
    desc(normalized_rank),
    desc(minus_log10_qvalue),
    desc(signal_value),
    chrom,
    summit,
    sample_id
  ) |>
  mutate(
    candidate_id = row_number(),
    candidate_start = summit - 250L,
    candidate_end = summit + 250L
  )

candidate_gr_final <- GenomicRanges::GRanges(
  seqnames = atac_peak_candidates_final$chrom,
  ranges = IRanges::IRanges(
    start = atac_peak_candidates_final$candidate_start,
    end = atac_peak_candidates_final$candidate_end
  )
)

keep_candidate_final <- rep(
  FALSE,
  length(candidate_gr_final)
)

available_final <- rep(
  TRUE,
  length(candidate_gr_final)
)

for (i in seq_along(candidate_gr_final)) {
  
  if (!available_final[i]) {
    next
  }
  
  keep_candidate_final[i] <- TRUE
  
  overlapping <- GenomicRanges::findOverlaps(
    candidate_gr_final[i],
    candidate_gr_final,
    ignore.strand = TRUE
  )
  
  overlapping_indices <-
    S4Vectors::subjectHits(overlapping)
  
  available_final[overlapping_indices] <- FALSE
}

atac_refined_loci <- atac_peak_candidates_final |>
  filter(keep_candidate_final) |>
  arrange(
    chrom,
    candidate_start,
    candidate_end
  ) |>
  transmute(
    locus_id = row_number(),
    chrom,
    start = candidate_start,
    end = candidate_end,
    summit,
    representative_sample = sample_id,
    representative_group = group,
    representative_qvalue =
      minus_log10_qvalue,
    representative_signal =
      signal_value,
    representative_normalized_rank =
      normalized_rank
  )

cat("\nFinal refined locus count:\n")
print(nrow(atac_refined_loci))

cat("\nFinal refined locus width QC:\n")

print(
  atac_refined_loci |>
    summarise(
      n_loci = n(),
      min_width = min(end - start + 1L),
      median_width = median(end - start + 1L),
      max_width = max(end - start + 1L)
    )
)

# Confirm that final loci do not overlap.

refined_gr <- GenomicRanges::GRanges(
  seqnames = atac_refined_loci$chrom,
  ranges = IRanges::IRanges(
    start = atac_refined_loci$start,
    end = atac_refined_loci$end
  )
)

refined_hits_all <- GenomicRanges::findOverlaps(
  refined_gr,
  refined_gr,
  ignore.strand = TRUE
)

refined_hits_nonself <- refined_hits_all[
  S4Vectors::queryHits(refined_hits_all) <
    S4Vectors::subjectHits(refined_hits_all)
]

cat("\nOverlaps among final refined loci:\n")
print(length(refined_hits_nonself))

# Compare final coordinates with the exploratory Step 16b set.

cat("\nComparison with exploratory Step 16b set:\n")

print(
  tibble(
    version = c(
      "Step16b_percent_rank",
      "Step16d_deterministic_rank"
    ),
    n_loci = c(
      nrow(atac_archr_style_loci),
      nrow(atac_refined_loci)
    )
  )
)

# ------------------------------------------------------------
# 17a. Finalize narrowPeak summit coordinate convention
# ------------------------------------------------------------

atac_peaks_window_coordcheck <- atac_peaks_window |>
  mutate(
    summit_bed0 =
      start + summit_offset,
    summit_gr1 =
      start + summit_offset + 1L
  )

cat("\nCheck existing summit convention:\n")

print(
  atac_peaks_window_coordcheck |>
    summarise(
      n_peaks = n(),
      existing_equals_bed0 =
        sum(summit == summit_bed0),
      existing_equals_gr1 =
        sum(summit == summit_gr1),
      min_existing_minus_bed0 =
        min(summit - summit_bed0),
      max_existing_minus_bed0 =
        max(summit - summit_bed0)
    )
)

cat("\nExample coordinate conversion:\n")

print(
  atac_peaks_window_coordcheck |>
    select(
      chrom,
      start,
      end,
      summit_offset,
      summit,
      summit_bed0,
      summit_gr1,
      sample_id,
      group
    ) |>
    slice_head(n = 10),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 17b. Final coordinate-corrected refined ATAC loci
# ------------------------------------------------------------

atac_refined_loci_final <- atac_refined_loci |>
  mutate(
    # Existing summit is BED 0-based.
    summit_bed0 = summit,
    
    # Equivalent 1-based genomic coordinate.
    summit_gr1 = summit_bed0 + 1L,
    
    # 501-bp interval in 1-based inclusive coordinates.
    start_gr1 = summit_gr1 - 250L,
    end_gr1 = summit_gr1 + 250L,
    
    # Equivalent BED interval:
    # 0-based start, end-exclusive.
    start_bed0 = start_gr1 - 1L,
    end_bed = end_gr1
  ) |>
  select(
    locus_id,
    chrom,
    start_gr1,
    end_gr1,
    start_bed0,
    end_bed,
    summit_gr1,
    summit_bed0,
    representative_sample,
    representative_group,
    representative_qvalue,
    representative_signal,
    representative_normalized_rank
  )

cat("\nFinal coordinate-corrected ATAC loci:\n")
print(nrow(atac_refined_loci_final))

cat("\nWidth and coordinate QC:\n")

print(
  atac_refined_loci_final |>
    summarise(
      n_loci = n(),
      
      min_gr_width =
        min(end_gr1 - start_gr1 + 1L),
      
      max_gr_width =
        max(end_gr1 - start_gr1 + 1L),
      
      min_bed_width =
        min(end_bed - start_bed0),
      
      max_bed_width =
        max(end_bed - start_bed0),
      
      all_summits_centered =
        all(
          summit_gr1 - start_gr1 == 250L &
            end_gr1 - summit_gr1 == 250L
        ),
      
      all_bed_gr_roundtrip =
        all(
          start_bed0 == start_gr1 - 1L &
            end_bed == end_gr1
        )
    )
)

cat("\nFirst 10 final loci:\n")

print(
  atac_refined_loci_final |>
    slice_head(n = 10),
  n = Inf,
  width = Inf
)

# Final non-overlap check using the corrected 1-based coordinates.

final_loci_gr <- GenomicRanges::GRanges(
  seqnames = atac_refined_loci_final$chrom,
  ranges = IRanges::IRanges(
    start = atac_refined_loci_final$start_gr1,
    end = atac_refined_loci_final$end_gr1
  )
)

final_hits_all <- GenomicRanges::findOverlaps(
  final_loci_gr,
  final_loci_gr,
  ignore.strand = TRUE
)

final_hits_nonself <- final_hits_all[
  S4Vectors::queryHits(final_hits_all) <
    S4Vectors::subjectHits(final_hits_all)
]

cat("\nOverlaps among coordinate-corrected final loci:\n")
print(length(final_hits_nonself))

# ------------------------------------------------------------
# 18a. Patient peak presence across final refined ATAC loci
# ------------------------------------------------------------

# Final loci are already represented as 1-based inclusive GRanges.
loci_gr <- GenomicRanges::GRanges(
  seqnames = atac_refined_loci_final$chrom,
  ranges = IRanges::IRanges(
    start = atac_refined_loci_final$start_gr1,
    end = atac_refined_loci_final$end_gr1
  )
)

# Convert original narrowPeak BED intervals to
# 1-based inclusive GRanges:
#
# BED start -> +1
# BED end   -> unchanged
peaks_gr <- GenomicRanges::GRanges(
  seqnames = atac_peaks_window$chrom,
  ranges = IRanges::IRanges(
    start = atac_peaks_window$start + 1L,
    end = atac_peaks_window$end
  )
)

peak_locus_hits <- GenomicRanges::findOverlaps(
  loci_gr,
  peaks_gr,
  ignore.strand = TRUE
)

peak_hit_table <- tibble(
  locus_id =
    atac_refined_loci_final$locus_id[
      S4Vectors::queryHits(peak_locus_hits)
    ],
  sample_id =
    atac_peaks_window$sample_id[
      S4Vectors::subjectHits(peak_locus_hits)
    ],
  group =
    atac_peaks_window$group[
      S4Vectors::subjectHits(peak_locus_hits)
    ]
) |>
  count(
    locus_id,
    sample_id,
    group,
    name = "n_overlapping_peaks"
  )

# Explicit 151 x 16 patient matrix.
atac_patient_lookup <- atac_peaks_window |>
  distinct(
    sample_id,
    group
  ) |>
  arrange(
    group,
    sample_id
  )

atac_locus_patient_peak <- tidyr::crossing(
  locus_id = atac_refined_loci_final$locus_id,
  atac_patient_lookup
) |>
  left_join(
    peak_hit_table,
    by = c(
      "locus_id",
      "sample_id",
      "group"
    )
  ) |>
  mutate(
    n_overlapping_peaks =
      tidyr::replace_na(
        n_overlapping_peaks,
        0L
      ),
    peak_present =
      n_overlapping_peaks > 0L
  )

cat("\nFinal locus x patient peak matrix dimensions:\n")
print(dim(atac_locus_patient_peak))

cat("\nPatients represented per locus:\n")
print(
  atac_locus_patient_peak |>
    count(locus_id) |>
    summarise(
      min_patients = min(n),
      max_patients = max(n)
    )
)

cat("\nPeak-presence totals by biological group:\n")
print(
  atac_locus_patient_peak |>
    group_by(
      locus_id,
      group
    ) |>
    summarise(
      n_patients = n(),
      n_peak_positive =
        sum(peak_present),
      peak_fraction =
        mean(peak_present),
      .groups = "drop"
    ) |>
    group_by(group) |>
    summarise(
      n_loci = n(),
      min_n_patients =
        min(n_patients),
      max_n_patients =
        max(n_patients),
      min_peak_fraction =
        min(peak_fraction),
      median_peak_fraction =
        median(peak_fraction),
      max_peak_fraction =
        max(peak_fraction),
      .groups = "drop"
    )
)

cat("\nDistribution of number of overlapping MACS2 peaks per patient/locus:\n")
print(
  atac_locus_patient_peak |>
    count(
      n_overlapping_peaks,
      name = "n_locus_patient_combinations"
    ) |>
    arrange(n_overlapping_peaks),
  n = Inf
)

# ------------------------------------------------------------
# 18b. Quantitative RPGC signal across final refined ATAC loci
# ------------------------------------------------------------

# atac_10bp uses BED-style genomic intervals:
# start = 0-based
# end   = end-exclusive.
#
# Convert the final loci to the same BED coordinate system.

final_loci_bed <- atac_refined_loci_final |>
  select(
    locus_id,
    chrom,
    locus_start = start_bed0,
    locus_end = end_bed
  )

# Restrict 10-bp signal bins to those that overlap at least
# one final locus.

signal_gr <- GenomicRanges::GRanges(
  seqnames = atac_10bp$chrom,
  ranges = IRanges::IRanges(
    start = atac_10bp$start + 1L,
    end = atac_10bp$end
  )
)

locus_gr_signal <- GenomicRanges::GRanges(
  seqnames = final_loci_bed$chrom,
  ranges = IRanges::IRanges(
    start = final_loci_bed$locus_start + 1L,
    end = final_loci_bed$locus_end
  )
)

signal_hits <- GenomicRanges::findOverlaps(
  locus_gr_signal,
  signal_gr,
  ignore.strand = TRUE
)

signal_query <- S4Vectors::queryHits(signal_hits)
signal_subject <- S4Vectors::subjectHits(signal_hits)

# Calculate exact overlap width in base pairs for each
# locus x signal-bin overlap.

overlap_start <- pmax(
  GenomicRanges::start(locus_gr_signal)[signal_query],
  GenomicRanges::start(signal_gr)[signal_subject]
)

overlap_end <- pmin(
  GenomicRanges::end(locus_gr_signal)[signal_query],
  GenomicRanges::end(signal_gr)[signal_subject]
)

overlap_width <- overlap_end - overlap_start + 1L

atac_locus_signal_hits <- tibble(
  locus_id =
    final_loci_bed$locus_id[signal_query],
  sample_id =
    atac_10bp$sample_id[signal_subject],
  group =
    atac_10bp$group[signal_subject],
  value =
    atac_10bp$value[signal_subject],
  overlap_width =
    overlap_width
)

# Weighted mean RPGC for each patient across each 501-bp locus.

atac_locus_patient_signal <- atac_locus_signal_hits |>
  group_by(
    locus_id,
    sample_id,
    group
  ) |>
  summarise(
    covered_bp =
      sum(overlap_width),
    mean_RPGC =
      weighted.mean(
        value,
        w = overlap_width
      ),
    .groups = "drop"
  )

cat("\nLocus x patient RPGC matrix dimensions:\n")
print(dim(atac_locus_patient_signal))

cat("\nCoverage QC across 501-bp loci:\n")

print(
  atac_locus_patient_signal |>
    summarise(
      n_rows = n(),
      min_covered_bp =
        min(covered_bp),
      median_covered_bp =
        median(covered_bp),
      max_covered_bp =
        max(covered_bp),
      min_mean_RPGC =
        min(mean_RPGC),
      median_mean_RPGC =
        median(mean_RPGC),
      max_mean_RPGC =
        max(mean_RPGC)
    )
)

# Join quantitative signal to the peak-presence matrix.

atac_locus_patient_final <- atac_locus_patient_peak |>
  left_join(
    atac_locus_patient_signal,
    by = c(
      "locus_id",
      "sample_id",
      "group"
    )
  )

cat("\nCombined peak + RPGC matrix QC:\n")

print(
  atac_locus_patient_final |>
    summarise(
      n_rows = n(),
      n_loci =
        n_distinct(locus_id),
      n_patients =
        n_distinct(sample_id),
      n_missing_signal =
        sum(is.na(mean_RPGC)),
      min_covered_bp =
        min(covered_bp, na.rm = TRUE),
      max_covered_bp =
        max(covered_bp, na.rm = TRUE)
    )
)

cat("\nPeak-positive versus peak-negative RPGC:\n")

print(
  atac_locus_patient_final |>
    group_by(peak_present) |>
    summarise(
      n = n(),
      median_RPGC =
        median(mean_RPGC),
      mean_RPGC =
        mean(mean_RPGC),
      p90_RPGC =
        quantile(
          mean_RPGC,
          0.90
        ),
      .groups = "drop"
    )
)

# ------------------------------------------------------------
# 19. AML/PBMC summaries across final refined ATAC loci
# ------------------------------------------------------------

atac_locus_group_final <- atac_locus_patient_final |>
  group_by(
    locus_id,
    group
  ) |>
  summarise(
    n_patients = n(),
    n_peak_positive =
      sum(peak_present),
    peak_fraction =
      mean(peak_present),
    
    group_mean_RPGC =
      mean(mean_RPGC),
    
    group_median_RPGC =
      median(mean_RPGC),
    
    .groups = "drop"
  )

atac_locus_group_wide <- atac_locus_group_final |>
  tidyr::pivot_wider(
    id_cols = locus_id,
    names_from = group,
    values_from = c(
      n_patients,
      n_peak_positive,
      peak_fraction,
      group_mean_RPGC,
      group_median_RPGC
    ),
    names_glue = "{group}_{.value}"
  )

atac_locus_summary_final <- atac_refined_loci_final |>
  left_join(
    atac_locus_group_wide,
    by = "locus_id"
  ) |>
  mutate(
    peak_fraction_difference =
      AML_peak_fraction -
      PBMC_peak_fraction,
    
    mean_RPGC_difference =
      AML_group_mean_RPGC -
      PBMC_group_mean_RPGC,
    
    median_RPGC_difference =
      AML_group_median_RPGC -
      PBMC_group_median_RPGC
  )

cat("\nCorrected group-summary QC:\n")

print(
  atac_locus_summary_final |>
    summarise(
      n_loci = n(),
      AML_n_min =
        min(AML_n_patients),
      AML_n_max =
        max(AML_n_patients),
      PBMC_n_min =
        min(PBMC_n_patients),
      PBMC_n_max =
        max(PBMC_n_patients),
      n_missing_mean =
        sum(
          is.na(AML_group_mean_RPGC) |
            is.na(PBMC_group_mean_RPGC)
        ),
      n_missing_median =
        sum(
          is.na(AML_group_median_RPGC) |
            is.na(PBMC_group_median_RPGC)
        )
    )
)

cat("\nMean versus median are now independently calculated:\n")

print(
  atac_locus_summary_final |>
    select(
      locus_id,
      AML_group_mean_RPGC,
      AML_group_median_RPGC,
      PBMC_group_mean_RPGC,
      PBMC_group_median_RPGC,
      mean_RPGC_difference,
      median_RPGC_difference
    ) |>
    slice_head(n = 10),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 20. Inspect final ATAC loci at DNMT3A and exact intron 2
# ------------------------------------------------------------

dnmt3a_start <- 25227874L
dnmt3a_end <- 25342590L

intron2_start <- 25300244L
intron2_end <- 25313912L

atac_dnmt3a_loci_final <- atac_locus_summary_final |>
  filter(
    chrom == "chr2",
    end_gr1 >= dnmt3a_start,
    start_gr1 <= dnmt3a_end
  ) |>
  mutate(
    overlaps_intron2 =
      end_gr1 >= intron2_start &
      start_gr1 <= intron2_end,
    
    summit_in_intron2 =
      summit_gr1 >= intron2_start &
      summit_gr1 <= intron2_end
  ) |>
  arrange(start_gr1)

atac_intron2_loci_final <- atac_dnmt3a_loci_final |>
  filter(overlaps_intron2)

cat("\nNumber of final loci overlapping DNMT3A:\n")
print(nrow(atac_dnmt3a_loci_final))

cat("\nNumber of final loci overlapping exact intron 2:\n")
print(nrow(atac_intron2_loci_final))

cat("\nFinal DNMT3A loci:\n")

print(
  atac_dnmt3a_loci_final |>
    select(
      locus_id,
      start_gr1,
      end_gr1,
      summit_gr1,
      overlaps_intron2,
      summit_in_intron2,
      AML_n_peak_positive,
      PBMC_n_peak_positive,
      AML_peak_fraction,
      PBMC_peak_fraction,
      AML_group_mean_RPGC,
      PBMC_group_mean_RPGC,
      mean_RPGC_difference,
      AML_group_median_RPGC,
      PBMC_group_median_RPGC,
      median_RPGC_difference
    ),
  n = Inf,
  width = Inf
)

cat("\nFinal exact-intron-2 loci:\n")

print(
  atac_intron2_loci_final |>
    select(
      locus_id,
      start_gr1,
      end_gr1,
      summit_gr1,
      representative_sample,
      representative_group,
      AML_n_peak_positive,
      PBMC_n_peak_positive,
      AML_peak_fraction,
      PBMC_peak_fraction,
      AML_group_mean_RPGC,
      PBMC_group_mean_RPGC,
      mean_RPGC_difference,
      AML_group_median_RPGC,
      PBMC_group_median_RPGC,
      median_RPGC_difference
    ),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 21. Save finalized ATAC locus metadata
# ------------------------------------------------------------

readr::write_tsv(
  atac_refined_loci_final,
  here::here(
    "metadata",
    "ATAC_DNMT3A_window_refined_loci.tsv"
  )
)

readr::write_tsv(
  atac_locus_patient_final,
  here::here(
    "metadata",
    "ATAC_DNMT3A_window_locus_patient_matrix.tsv"
  )
)

readr::write_tsv(
  atac_locus_summary_final,
  here::here(
    "metadata",
    "ATAC_DNMT3A_window_locus_group_summary.tsv"
  )
)

readr::write_tsv(
  atac_dnmt3a_loci_final,
  here::here(
    "metadata",
    "ATAC_DNMT3A_gene_loci.tsv"
  )
)

readr::write_tsv(
  atac_intron2_loci_final,
  here::here(
    "metadata",
    "ATAC_DNMT3A_intron2_loci.tsv"
  )
)

cat("\nSaved finalized ATAC metadata files:\n")

atac_saved_files <- c(
  "ATAC_DNMT3A_window_refined_loci.tsv",
  "ATAC_DNMT3A_window_locus_patient_matrix.tsv",
  "ATAC_DNMT3A_window_locus_group_summary.tsv",
  "ATAC_DNMT3A_gene_loci.tsv",
  "ATAC_DNMT3A_intron2_loci.tsv"
)

atac_saved_paths <- here::here(
  "metadata",
  atac_saved_files
)

print(
  tibble(
    file = atac_saved_files,
    exists = file.exists(atac_saved_paths),
    size_bytes = file.info(atac_saved_paths)$size
  ),
  n = Inf
)

cat("\nFinal saved-row QC:\n")

print(
  tibble(
    dataset = c(
      "refined loci",
      "locus x patient",
      "locus group summary",
      "DNMT3A loci",
      "intron 2 loci"
    ),
    expected_rows = c(
      151L,
      2416L,
      151L,
      23L,
      2L
    ),
    saved_rows = c(
      nrow(
        readr::read_tsv(
          atac_saved_paths[1],
          show_col_types = FALSE
        )
      ),
      nrow(
        readr::read_tsv(
          atac_saved_paths[2],
          show_col_types = FALSE
        )
      ),
      nrow(
        readr::read_tsv(
          atac_saved_paths[3],
          show_col_types = FALSE
        )
      ),
      nrow(
        readr::read_tsv(
          atac_saved_paths[4],
          show_col_types = FALSE
        )
      ),
      nrow(
        readr::read_tsv(
          atac_saved_paths[5],
          show_col_types = FALSE
        )
      )
    )
  )
)

# ------------------------------------------------------------
# 22. Build AML and PBMC mean ATAC RPGC tracks
# ------------------------------------------------------------

atac_group_mean_10bp <- atac_10bp |>
  group_by(
    chrom,
    start,
    end,
    group
  ) |>
  summarise(
    n_patients = n(),
    mean_RPGC = mean(value),
    .groups = "drop"
  )

cat("\nGroup mean track dimensions:\n")
print(dim(atac_group_mean_10bp))

cat("\nGroup mean track QC:\n")

print(
  atac_group_mean_10bp |>
    group_by(group) |>
    summarise(
      n_bins = n(),
      min_n_patients = min(n_patients),
      max_n_patients = max(n_patients),
      min_RPGC = min(mean_RPGC),
      median_RPGC = median(mean_RPGC),
      mean_RPGC_over_window = mean(mean_RPGC),
      p95_RPGC = quantile(mean_RPGC, 0.95),
      p99_RPGC = quantile(mean_RPGC, 0.99),
      p995_RPGC = quantile(mean_RPGC, 0.995),
      p999_RPGC = quantile(mean_RPGC, 0.999),
      max_RPGC = max(mean_RPGC),
      .groups = "drop"
    )
)

cat("\nExpected patient counts:\n")

print(
  atac_group_mean_10bp |>
    count(
      group,
      n_patients
    )
)

cat("\nCoordinate-grid QC:\n")

print(
  atac_group_mean_10bp |>
    group_by(group) |>
    summarise(
      first_start = min(start),
      last_end = max(end),
      min_width = min(end - start),
      max_width = max(end - start),
      n_bins = n(),
      .groups = "drop"
    )
)

cat("\nAML/PBMC coordinate grids identical:\n")

atac_aml_grid <- atac_group_mean_10bp |>
  filter(group == "AML") |>
  select(chrom, start, end)

atac_pbmc_grid <- atac_group_mean_10bp |>
  filter(group == "PBMC") |>
  select(chrom, start, end)

print(
  identical(
    atac_aml_grid,
    atac_pbmc_grid
  )
)

cat("\nLargest group-mean signal values:\n")

print(
  atac_group_mean_10bp |>
    arrange(desc(mean_RPGC)) |>
    select(
      chrom,
      start,
      end,
      group,
      mean_RPGC
    ) |>
    slice_head(n = 20),
  n = Inf
)

# ------------------------------------------------------------
# 23. Save AML and PBMC mean ATAC tracks as bedGraph
# ------------------------------------------------------------

atac_group_track_dir <- here::here(
  "tracks",
  "ATAC_group_summary"
)

dir.create(
  atac_group_track_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write_atac_group_bedgraph <- function(group_name) {
  
  output_file <- here::here(
    "tracks",
    "ATAC_group_summary",
    paste0(
      "ATAC_",
      group_name,
      "_mean_RPGC_DNMT3A.bedGraph"
    )
  )
  
  track_data <- atac_group_mean_10bp |>
    filter(group == group_name) |>
    arrange(chrom, start) |>
    select(
      chrom,
      start,
      end,
      mean_RPGC
    )
  
  readr::write_tsv(
    track_data,
    output_file,
    col_names = FALSE
  )
  
  tibble(
    group = group_name,
    file = output_file,
    n_rows = nrow(track_data),
    min_value = min(track_data$mean_RPGC),
    max_value = max(track_data$mean_RPGC),
    size_bytes = file.info(output_file)$size
  )
}

atac_group_bedgraph_status <- bind_rows(
  write_atac_group_bedgraph("AML"),
  write_atac_group_bedgraph("PBMC")
)

cat("\nATAC group bedGraph files:\n")

print(
  atac_group_bedgraph_status,
  n = Inf,
  width = Inf
)

cat("\nFirst five AML bedGraph rows:\n")

print(
  readr::read_tsv(
    atac_group_bedgraph_status$file[
      atac_group_bedgraph_status$group == "AML"
    ],
    col_names = c(
      "chrom",
      "start",
      "end",
      "value"
    ),
    show_col_types = FALSE
  ) |>
    slice_head(n = 5)
)

cat("\nFirst five PBMC bedGraph rows:\n")

print(
  readr::read_tsv(
    atac_group_bedgraph_status$file[
      atac_group_bedgraph_status$group == "PBMC"
    ],
    col_names = c(
      "chrom",
      "start",
      "end",
      "value"
    ),
    show_col_types = FALSE
  ) |>
    slice_head(n = 5)
)

# ------------------------------------------------------------
# 24. Convert AML and PBMC mean ATAC bedGraphs to BigWig
# ------------------------------------------------------------

windows_to_wsl <- function(path) {
  
  path <- normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )
  
  drive <- tolower(substr(path, 1, 1))
  rest <- substr(path, 3, nchar(path))
  
  paste0(
    "/mnt/",
    drive,
    rest
  )
}

atac_group_bigwig_status <- atac_group_bedgraph_status |>
  rowwise() |>
  mutate(
    bigwig_file =
      sub(
        "\\.bedGraph$",
        ".bw",
        file
      ),
    
    bedgraph_wsl =
      windows_to_wsl(file),
    
    bigwig_wsl =
      windows_to_wsl(bigwig_file),
    
    command =
      paste(
        "~/ucsc_tools/bedGraphToBigWig",
        shQuote(bedgraph_wsl),
        "~/ucsc_tools/hg38.chrom.sizes",
        shQuote(bigwig_wsl)
      ),
    
    exit_status =
      system2(
        "wsl",
        c(
          "bash",
          "-lc",
          shQuote(command)
        )
      )
  ) |>
  ungroup() |>
  mutate(
    bigwig_exists =
      file.exists(bigwig_file),
    
    bigwig_size_bytes =
      ifelse(
        bigwig_exists,
        file.info(bigwig_file)$size,
        NA_real_
      )
  )

cat("\nATAC group BigWig conversion status:\n")

print(
  atac_group_bigwig_status |>
    select(
      group,
      exit_status,
      bigwig_exists,
      bigwig_size_bytes,
      bigwig_file
    ),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 24b. Validate final ATAC group BigWigs
# ------------------------------------------------------------

validate_atac_group_bigwig <- function(group_name) {
  
  bw_file <- atac_group_bigwig_status |>
    filter(group == group_name) |>
    pull(bigwig_file)
  
  bw <- cpp11bigwig::read_bigwig(
    bw_file,
    chrom = "chr2",
    start = 24727880,
    end = 25842590
  )
  
  tibble(
    group = group_name,
    n_intervals = nrow(bw),
    first_start = min(bw$start),
    last_end = max(bw$end),
    min_value = min(bw$value),
    max_value = max(bw$value)
  )
}

atac_group_bigwig_qc <- bind_rows(
  validate_atac_group_bigwig("AML"),
  validate_atac_group_bigwig("PBMC")
)

cat("\nATAC group BigWig validation:\n")

print(
  atac_group_bigwig_qc,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 25. Prepare derived ATAC loci for UCSC bigBed
# ------------------------------------------------------------

atac_loci_bed <- atac_locus_summary_final |>
  transmute(
    chrom = chrom,
    chromStart = start_bed0,
    chromEnd = end_bed,
    
    name =
      paste0(
        "ATAC_locus_",
        locus_id
      ),
    
    score =
      round(
        1000 *
          pmax(
            AML_peak_fraction,
            PBMC_peak_fraction
          )
      ),
    
    strand = ".",
    
    thickStart = summit_bed0,
    thickEnd = summit_bed0 + 1L,
    
    itemRgb = "80,80,80",
    
    locus_id = locus_id,
    
    AML_peak_n =
      AML_n_peak_positive,
    
    AML_total_n =
      AML_n_patients,
    
    AML_peak_fraction =
      AML_peak_fraction,
    
    PBMC_peak_n =
      PBMC_n_peak_positive,
    
    PBMC_total_n =
      PBMC_n_patients,
    
    PBMC_peak_fraction =
      PBMC_peak_fraction,
    
    AML_mean_RPGC =
      AML_group_mean_RPGC,
    
    PBMC_mean_RPGC =
      PBMC_group_mean_RPGC
  ) |>
  arrange(
    chrom,
    chromStart,
    chromEnd
  )

cat("\nDerived ATAC UCSC locus table dimensions:\n")
print(dim(atac_loci_bed))

cat("\nBED coordinate and score QC:\n")

print(
  atac_loci_bed |>
    summarise(
      n_loci = n(),
      min_width =
        min(chromEnd - chromStart),
      max_width =
        max(chromEnd - chromStart),
      min_score =
        min(score),
      max_score =
        max(score),
      min_AML_fraction =
        min(AML_peak_fraction),
      max_AML_fraction =
        max(AML_peak_fraction),
      min_PBMC_fraction =
        min(PBMC_peak_fraction),
      max_PBMC_fraction =
        max(PBMC_peak_fraction)
    )
)

cat("\nFirst 10 derived ATAC UCSC loci:\n")

print(
  atac_loci_bed |>
    slice_head(n = 10),
  n = Inf,
  width = Inf
)

cat("\nExact intron-2 derived loci:\n")

print(
  atac_loci_bed |>
    filter(
      chromEnd > 25300243L,
      chromStart < 25313912L
    ),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 26. Write final derived ATAC locus BED9+9
# ------------------------------------------------------------

atac_loci_bed_final <- atac_loci_bed |>
  mutate(
    score = 0L
  )

atac_loci_bed_file <- here::here(
  "tracks",
  "ATAC_group_summary",
  "ATAC_derived_501bp_loci_DNMT3A.bed"
)

readr::write_tsv(
  atac_loci_bed_final,
  atac_loci_bed_file,
  col_names = FALSE
)

cat("\nFinal derived ATAC BED file:\n")
print(atac_loci_bed_file)

cat("\nFinal BED file QC:\n")

print(
  tibble(
    exists = file.exists(atac_loci_bed_file),
    size_bytes = file.info(atac_loci_bed_file)$size,
    n_rows = nrow(atac_loci_bed_final),
    n_columns = ncol(atac_loci_bed_final),
    min_score = min(atac_loci_bed_final$score),
    max_score = max(atac_loci_bed_final$score)
  )
)

cat("\nFirst three final BED rows:\n")

print(
  atac_loci_bed_final |>
    slice_head(n = 3),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 27. Create AutoSql schema for derived ATAC loci
# ------------------------------------------------------------

atac_loci_as_file <- here::here(
  "tracks",
  "ATAC_group_summary",
  "ATAC_derived_501bp_loci_DNMT3A.as"
)

atac_loci_as <- c(
  "table ATACDerivedLoci",
  "\"Derived fixed-width ATAC loci with AML/PBMC descriptive summaries\"",
  "(",
  "    string chrom;              \"Reference sequence chromosome\"",
  "    uint chromStart;           \"Start position in chromosome (0-based)\"",
  "    uint chromEnd;             \"End position in chromosome\"",
  "    string name;               \"Derived ATAC locus identifier\"",
  "    uint score;                \"Neutral BED score\"",
  "    char[1] strand;            \"Strand; not applicable for ATAC loci\"",
  "    uint thickStart;           \"Representative summit position\"",
  "    uint thickEnd;             \"Representative summit position plus one\"",
  "    uint reserved;             \"Item RGB value\"",
  "    uint locus_id;             \"Derived locus numeric identifier\"",
  "    uint AML_peak_n;           \"AML patients with an overlapping deposited MACS2 peak\"",
  "    uint AML_total_n;          \"Total AML patients\"",
  "    float AML_peak_fraction;   \"Fraction of AML patients with an overlapping deposited MACS2 peak\"",
  "    uint PBMC_peak_n;          \"PBMC samples with an overlapping deposited MACS2 peak\"",
  "    uint PBMC_total_n;         \"Total PBMC samples\"",
  "    float PBMC_peak_fraction;  \"Fraction of PBMC samples with an overlapping deposited MACS2 peak\"",
  "    float AML_mean_RPGC;       \"Mean ATAC RPGC signal across AML patients over the 501-bp locus\"",
  "    float PBMC_mean_RPGC;      \"Mean ATAC RPGC signal across PBMC samples over the 501-bp locus\"",
  ")"
)

writeLines(
  atac_loci_as,
  atac_loci_as_file
)

cat("\nAutoSql schema file:\n")
print(atac_loci_as_file)

cat("\nAutoSql schema contents:\n")
cat(
  readLines(atac_loci_as_file),
  sep = "\n"
)

cat("\n\nSchema QC:\n")

print(
  tibble(
    exists =
      file.exists(atac_loci_as_file),
    size_bytes =
      file.info(atac_loci_as_file)$size,
    n_schema_fields =
      sum(
        grepl(
          ";",
          readLines(atac_loci_as_file),
          fixed = TRUE
        )
      ),
    n_bed_columns =
      ncol(atac_loci_bed_final)
  )
)

# ------------------------------------------------------------
# 28. Convert derived ATAC loci to annotated bigBed
# ------------------------------------------------------------

atac_loci_bed_bigbed <- atac_loci_bed_final |>
  mutate(
    itemRgb = 5263440L
  )

readr::write_tsv(
  atac_loci_bed_bigbed,
  atac_loci_bed_file,
  col_names = FALSE
)

atac_loci_bb_file <- sub(
  "\\.bed$",
  ".bb",
  atac_loci_bed_file
)

atac_loci_bed_wsl <-
  windows_to_wsl(atac_loci_bed_file)

atac_loci_as_wsl <-
  windows_to_wsl(atac_loci_as_file)

atac_loci_bb_wsl <-
  windows_to_wsl(atac_loci_bb_file)

atac_bigbed_command <- paste0(
  "~/ucsc_tools/bedToBigBed ",
  "-as=", shQuote(atac_loci_as_wsl), " ",
  "-type=bed9+9 ",
  "-tab ",
  shQuote(atac_loci_bed_wsl), " ",
  "~/ucsc_tools/hg38.chrom.sizes ",
  shQuote(atac_loci_bb_wsl)
)

cat("\nRunning bedToBigBed:\n")
cat(atac_bigbed_command, "\n")

atac_bigbed_exit_status <- system2(
  "wsl",
  c(
    "bash",
    "-lc",
    shQuote(atac_bigbed_command)
  )
)

cat("\nDerived ATAC bigBed conversion QC:\n")

print(
  tibble(
    exit_status =
      atac_bigbed_exit_status,
    bigbed_exists =
      file.exists(atac_loci_bb_file),
    size_bytes =
      ifelse(
        file.exists(atac_loci_bb_file),
        file.info(atac_loci_bb_file)$size,
        NA_real_
      ),
    file =
      atac_loci_bb_file
  )
)

# ------------------------------------------------------------
# 28b. Validate derived ATAC bigBed with UCSC bigBedInfo
# ------------------------------------------------------------

atac_bigbed_info_command <- paste(
  "~/ucsc_tools/bigBedInfo",
  shQuote(atac_loci_bb_wsl)
)

atac_bigbed_info <- system2(
  "wsl",
  c(
    "bash",
    "-lc",
    shQuote(atac_bigbed_info_command)
  ),
  stdout = TRUE,
  stderr = TRUE
)

cat("\nUCSC bigBedInfo output:\n")

cat(
  atac_bigbed_info,
  sep = "\n"
)