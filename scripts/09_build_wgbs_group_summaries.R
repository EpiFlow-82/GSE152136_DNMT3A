library(tidyverse)
library(here)
library(cpp11bigwig)

# ============================================================
# WGBS group summaries
# Step 1: validate 10-bp coordinate structure
# ============================================================

chrom <- "chr2"

# DNMT3A +/- 500 kb, 1-based inclusive
window_start <- 24727874
window_end   <- 25842590

# ============================================================
# Load WGBS manifest
# ============================================================

bigwig_manifest <- read_tsv(
  here("metadata", "regional_1d_bigwig_tracks.tsv"),
  show_col_types = FALSE
)

wgbs_tracks <- bigwig_manifest |>
  filter(track_type == "WGBS") |>
  arrange(group, sample_id)

cat("WGBS tracks:", nrow(wgbs_tracks), "\n")

wgbs_tracks |>
  count(group) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ============================================================
# Inspect interval coordinates for all 20 WGBS tracks
# ============================================================

grid_check <- lapply(
  seq_len(nrow(wgbs_tracks)),
  function(i) {
    
    x <- wgbs_tracks[i, ]
    
    bw <- cpp11bigwig::read_bigwig(
      x$bigwig_file,
      chrom = chrom,
      start = window_start,
      end = window_end
    )
    
    tibble(
      sample_id = x$sample_id,
      group = x$group,
      n_intervals = nrow(bw),
      
      # BigWig coordinates returned by cpp11bigwig are
      # 0-based start, end-exclusive.
      all_widths_multiple_10 =
        all((bw$end - bw$start) %% 10 == 0),
      
      start_mod10_values =
        paste(sort(unique(bw$start %% 10)), collapse = ","),
      
      end_mod10_values =
        paste(sort(unique(bw$end %% 10)), collapse = ","),
      
      min_start = min(bw$start),
      max_end = max(bw$end)
    )
  }
) |>
  bind_rows()

cat("\n10-bp grid check for all WGBS samples:\n")

grid_check |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nUnique start modulo-10 patterns:\n")

grid_check |>
  count(start_mod10_values) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nUnique end modulo-10 patterns:\n")

grid_check |>
  count(end_mod10_values) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nAll interval widths multiples of 10 bp:",
    all(grid_check$all_widths_multiple_10),
    "\n")

# ============================================================
# Step 2: reconstruct every WGBS sample on common 10-bp grid
# ============================================================

grid_start <- min(grid_check$min_start)
grid_end   <- max(grid_check$max_end)

wgbs_grid <- tibble(
  chrom = chrom,
  start = seq(
    from = grid_start,
    to = grid_end - 10,
    by = 10
  )
) |>
  mutate(
    end = start + 10,
    bin_id = row_number()
  )

cat("\nCommon WGBS grid:\n")
cat("Start:", grid_start, "\n")
cat("End:", grid_end, "\n")
cat("10-bp bins:", nrow(wgbs_grid), "\n")

# ============================================================
# Function to expand one stored BigWig interval into 10-bp bins
# ============================================================

expand_wgbs_sample <- function(i) {
  
  x <- wgbs_tracks[i, ]
  
  bw <- cpp11bigwig::read_bigwig(
    x$bigwig_file,
    chrom = chrom,
    start = window_start,
    end = window_end
  )
  
  expanded <- lapply(
    seq_len(nrow(bw)),
    function(j) {
      
      starts <- seq(
        from = bw$start[j],
        to = bw$end[j] - 10,
        by = 10
      )
      
      tibble(
        chrom = chrom,
        start = starts,
        end = starts + 10,
        value = bw$value[j]
      )
    }
  ) |>
    bind_rows() |>
    mutate(
      sample_id = x$sample_id,
      group = x$group
    )
  
  expanded
}

# ============================================================
# Expand all 20 WGBS samples
# ============================================================

wgbs_10bp <- lapply(
  seq_len(nrow(wgbs_tracks)),
  expand_wgbs_sample
) |>
  bind_rows()

cat("\nExpanded 10-bp observations:",
    nrow(wgbs_10bp), "\n")

cat("\nExpanded bins per sample:\n")

expanded_counts <- wgbs_10bp |>
  count(group, sample_id, name = "n_observed_bins") |>
  arrange(group, sample_id)

expanded_counts |>
  as.data.frame() |>
  print(row.names = FALSE)

# ============================================================
# Validation
# ============================================================

cat("\nAll reconstructed intervals exactly 10 bp:",
    all(wgbs_10bp$end - wgbs_10bp$start == 10),
    "\n")

cat(
  "Duplicate sample/bin combinations:",
  wgbs_10bp |>
    count(sample_id, start, end) |>
    filter(n > 1) |>
    nrow(),
  "\n"
)

cat(
  "Signal range after reconstruction:",
  min(wgbs_10bp$value, na.rm = TRUE),
  "to",
  max(wgbs_10bp$value, na.rm = TRUE),
  "\n"
)

# Save reconstructed sample-level 10-bp data
write_tsv(
  wgbs_10bp,
  here(
    "metadata",
    "WGBS_sample_level_10bp.tsv"
  )
)

cat(
  "\nSaved sample-level 10-bp WGBS table:",
  here("metadata", "WGBS_sample_level_10bp.tsv"),
  "\n"
)
# ============================================================
# Step 3: calculate WGBS group summaries
#
# IMPORTANT:
# - Only observed values contribute to the mean.
# - Missing bins are NOT converted to zero.
# - n_samples records how many samples contributed.
# ============================================================

wgbs_group_summary <- wgbs_10bp |>
  group_by(group, chrom, start, end) |>
  summarise(
    mean_methylation = mean(value),
    n_samples = n(),
    .groups = "drop"
  ) |>
  arrange(group, chrom, start)

cat("\nGroup-summary bins:\n")

wgbs_group_summary |>
  count(group, name = "n_bins") |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nNumber of contributing samples per bin:\n")

wgbs_group_summary |>
  count(group, n_samples, name = "n_bins") |>
  arrange(group, n_samples) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nGroup mean methylation ranges:\n")

wgbs_group_summary |>
  summarise(
    minimum = min(mean_methylation),
    maximum = max(mean_methylation),
    .by = group
  ) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ============================================================
# Save detailed summary table
# ============================================================

write_tsv(
  wgbs_group_summary,
  here("metadata", "WGBS_group_summary_10bp.tsv")
)

cat("\nSaved:\n")
cat("metadata/WGBS_group_summary_10bp.tsv\n")

# ============================================================
# Step 4: apply minimum support for displayed group means
#
# AML:
#   require >= 16 of 18 samples
#
# PBMC:
#   require 2 of 2 samples
#
# The complete unfiltered summary remains preserved in:
# metadata/WGBS_group_summary_10bp.tsv
# ============================================================

wgbs_group_display <- wgbs_group_summary |>
  filter(
    (group == "AML"  & n_samples >= 16) |
      (group == "PBMC" & n_samples == 2)
  )

cat("\nDisplay-ready WGBS group summaries:\n")

display_counts <- wgbs_group_display |>
  count(group, name = "n_bins") |>
  mutate(
    nominal_group_n = if_else(group == "AML", 18L, 2L),
    minimum_required_n = if_else(group == "AML", 16L, 2L)
  )

display_counts |>
  as.data.frame() |>
  print(row.names = FALSE)

# ============================================================
# Save display-ready table
# ============================================================

write_tsv(
  wgbs_group_display,
  here("metadata", "WGBS_group_summary_10bp_display.tsv")
)

# ============================================================
# Create bedGraph files
#
# Coordinates are already:
# start = 0-based
# end   = end-exclusive
# ============================================================

dir.create(
  here("tracks", "group_summaries"),
  recursive = TRUE,
  showWarnings = FALSE
)

aml_bedgraph <- wgbs_group_display |>
  filter(group == "AML") |>
  select(
    chrom,
    start,
    end,
    mean_methylation
  )

pbmc_bedgraph <- wgbs_group_display |>
  filter(group == "PBMC") |>
  select(
    chrom,
    start,
    end,
    mean_methylation
  )

write_tsv(
  aml_bedgraph,
  here(
    "tracks",
    "group_summaries",
    "WGBS_AML_mean_n16of18_DNMT3A.bedGraph"
  ),
  col_names = FALSE
)

write_tsv(
  pbmc_bedgraph,
  here(
    "tracks",
    "group_summaries",
    "WGBS_PBMC_mean_n2of2_DNMT3A.bedGraph"
  ),
  col_names = FALSE
)

cat("\nSaved display table and bedGraphs:\n")

cat(
  "metadata/WGBS_group_summary_10bp_display.tsv\n"
)

cat(
  "tracks/group_summaries/WGBS_AML_mean_n16of18_DNMT3A.bedGraph\n"
)

cat(
  "tracks/group_summaries/WGBS_PBMC_mean_n2of2_DNMT3A.bedGraph\n"
)

# ============================================================
# Step 5: validate display-ready bedGraphs
# ============================================================

validate_bedgraph <- function(x, group_name) {
  
  cat("\n---", group_name, "---\n")
  
  cat("Rows:", nrow(x), "\n")
  
  cat(
    "All intervals exactly 10 bp:",
    all(x$end - x$start == 10),
    "\n"
  )
  
  cat(
    "Sorted by genomic start:",
    identical(x$start, sort(x$start)),
    "\n"
  )
  
  cat(
    "Overlapping intervals:",
    sum(x$start[-1] < x$end[-nrow(x)]),
    "\n"
  )
  
  cat(
    "Missing methylation values:",
    sum(is.na(x$mean_methylation)),
    "\n"
  )
  
  cat(
    "Values outside 0-1:",
    sum(
      x$mean_methylation < 0 |
        x$mean_methylation > 1,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "Methylation range:",
    min(x$mean_methylation),
    "to",
    max(x$mean_methylation),
    "\n"
  )
}

validate_bedgraph(
  aml_bedgraph,
  "AML mean >=16/18"
)

validate_bedgraph(
  pbmc_bedgraph,
  "PBMC mean 2/2"
)

# ============================================================
# Step 6: convert WGBS group-summary bedGraphs to BigWig
# using UCSC bedGraphToBigWig through WSL
# ============================================================

summary_dir <- here("tracks", "group_summaries")

aml_bg <- file.path(
  summary_dir,
  "WGBS_AML_mean_n16of18_DNMT3A.bedGraph"
)

pbmc_bg <- file.path(
  summary_dir,
  "WGBS_PBMC_mean_n2of2_DNMT3A.bedGraph"
)

aml_bw <- file.path(
  summary_dir,
  "WGBS_AML_mean_n16of18_DNMT3A.bw"
)

pbmc_bw <- file.path(
  summary_dir,
  "WGBS_PBMC_mean_n2of2_DNMT3A.bw"
)

# Convert Windows paths to WSL paths
wsl_path <- function(path) {
  result <- system2(
    "wsl",
    c("wslpath", "-a", shQuote(normalizePath(path))),
    stdout = TRUE
  )
  trimws(result)
}

aml_bg_wsl  <- wsl_path(aml_bg)
pbmc_bg_wsl <- wsl_path(pbmc_bg)

aml_bw_wsl <- wsl_path(
  normalizePath(
    dirname(aml_bw),
    mustWork = TRUE
  )
)
aml_bw_wsl <- paste0(
  aml_bw_wsl,
  "/",
  basename(aml_bw)
)

pbmc_bw_wsl <- wsl_path(
  normalizePath(
    dirname(pbmc_bw),
    mustWork = TRUE
  )
)
pbmc_bw_wsl <- paste0(
  pbmc_bw_wsl,
  "/",
  basename(pbmc_bw)
)

# UCSC resources already installed in Ubuntu/WSL
converter <- "~/ucsc_tools/bedGraphToBigWig"
chrom_sizes <- "~/ucsc_tools/hg38.chrom.sizes"

# Convert AML
cmd_aml <- paste(
  converter,
  shQuote(aml_bg_wsl),
  chrom_sizes,
  shQuote(aml_bw_wsl)
)

status_aml <- system2(
  "wsl",
  c("bash", "-lc", shQuote(cmd_aml))
)

# Convert PBMC
cmd_pbmc <- paste(
  converter,
  shQuote(pbmc_bg_wsl),
  chrom_sizes,
  shQuote(pbmc_bw_wsl)
)

status_pbmc <- system2(
  "wsl",
  c("bash", "-lc", shQuote(cmd_pbmc))
)

cat("\nConversion status:\n")
cat("AML :", status_aml, "\n")
cat("PBMC:", status_pbmc, "\n")

cat("\nBigWig files exist:\n")
cat("AML :", file.exists(aml_bw), "\n")
cat("PBMC:", file.exists(pbmc_bw), "\n")

cat("\nBigWig sizes:\n")
print(
  file.info(c(aml_bw, pbmc_bw))["size"]
)

# ============================================================
# Step 7: read back and validate generated BigWigs
# ============================================================

aml_bw_check <- cpp11bigwig::read_bigwig(
  aml_bw,
  chrom = chrom,
  start = window_start,
  end = window_end
)

pbmc_bw_check <- cpp11bigwig::read_bigwig(
  pbmc_bw,
  chrom = chrom,
  start = window_start,
  end = window_end
)

cat("\nBigWig read-back validation:\n")

readback_summary <- tibble(
  group = c("AML", "PBMC"),
  n_intervals = c(
    nrow(aml_bw_check),
    nrow(pbmc_bw_check)
  ),
  minimum = c(
    min(aml_bw_check$value),
    min(pbmc_bw_check$value)
  ),
  maximum = c(
    max(aml_bw_check$value),
    max(pbmc_bw_check$value)
  ),
  min_start = c(
    min(aml_bw_check$start),
    min(pbmc_bw_check$start)
  ),
  max_end = c(
    max(aml_bw_check$end),
    max(pbmc_bw_check$end)
  )
)

readback_summary |>
  as.data.frame() |>
  print(row.names = FALSE)

cat(
  "\nAML BigWig readable:",
  nrow(aml_bw_check) > 0,
  "\n"
)

cat(
  "PBMC BigWig readable:",
  nrow(pbmc_bw_check) > 0,
  "\n"
)

cat(
  "All read-back values within 0-1:",
  all(
    aml_bw_check$value >= 0 &
      aml_bw_check$value <= 1
  ) &&
    all(
      pbmc_bw_check$value >= 0 &
        pbmc_bw_check$value <= 1
    ),
  "\n"
)