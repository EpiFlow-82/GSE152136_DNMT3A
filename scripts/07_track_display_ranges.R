library(tidyverse)
library(here)
library(cpp11bigwig)

bigwig_manifest <- read_tsv(
  here("metadata", "regional_1d_bigwig_tracks.tsv"),
  show_col_types = FALSE
)

atac_tracks <- bigwig_manifest |>
  filter(track_type == "ATAC-seq")

atac_ranges <- lapply(
  seq_len(nrow(atac_tracks)),
  function(i) {
    
    x <- atac_tracks[i, ]
    
    bw <- cpp11bigwig::read_bigwig(
      x$bigwig_file,
      chrom = "chr2",
      start = 24727873,
      end = 25842590
    )
    
    tibble(
      sample_id = x$sample_id,
      group = x$group,
      n_intervals = nrow(bw),
      maximum = max(bw$value, na.rm = TRUE),
      p99 = quantile(bw$value, 0.99, na.rm = TRUE),
      p995 = quantile(bw$value, 0.995, na.rm = TRUE),
      p999 = quantile(bw$value, 0.999, na.rm = TRUE)
    )
  }
) |>
  bind_rows()

print(
  atac_ranges |>
    arrange(group, sample_id),
  n = Inf
)


wgbsm_tracks <- bigwig_manifest |>
  filter(track_type == "WGBS")

wgbsm_ranges <- lapply(
  seq_len(nrow(wgbsm_tracks)),
  function(i) {
    
    x <- wgbsm_tracks[i, ]
    
    bw <- cpp11bigwig::read_bigwig(
      x$bigwig_file,
      chrom = "chr2",
      start = 24727873,
      end = 25842590
    )
    
    tibble(
      sample_id = x$sample_id,
      group = x$group,
      n_intervals = nrow(bw),
      minimum = min(bw$value, na.rm = TRUE),
      maximum = max(bw$value, na.rm = TRUE)
    )
  }
) |>
  bind_rows()

print(
  wgbsm_ranges |>
    arrange(group, sample_id),
  n = Inf
)

cat("\nOverall WGBS value range:\n")
print(
  range(
    c(wgbsm_ranges$minimum, wgbsm_ranges$maximum),
    na.rm = TRUE
  )
)