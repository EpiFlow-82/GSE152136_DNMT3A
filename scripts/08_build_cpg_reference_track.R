library(tidyverse)
library(here)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(cpp11bigwig)
# ============================================================
# Coordinates
# ============================================================

chrom <- "chr2"

# DNMT3A +/- 500 kb window
window_start <- 24727874
window_end   <- 25842590

# DNMT3A intron 2, 1-based inclusive
intron2_start <- 25300244
intron2_end   <- 25313912

# ============================================================
# Extract hg38 reference sequence
# ============================================================

hg38 <- BSgenome.Hsapiens.UCSC.hg38

seq_window <- getSeq(
  hg38,
  names = chrom,
  start = window_start,
  end = window_end
)

cat("Window length:", length(seq_window), "bp\n")

# ============================================================
# Find all CpG dinucleotides
# ============================================================

cpg_hits <- matchPattern("CG", seq_window)

cpg <- tibble(
  chrom = chrom,
  start_1based = window_start + start(cpg_hits) - 1,
  end_1based   = window_start + end(cpg_hits) - 1
)

cat("CpGs in DNMT3A window:", nrow(cpg), "\n")

# ============================================================
# Identify CpGs in DNMT3A intron 2
# ============================================================

cpg <- cpg |>
  mutate(
    in_intron2 =
      start_1based >= intron2_start &
      end_1based <= intron2_end
  )

cpg_intron2 <- cpg |>
  filter(in_intron2)

cat("CpGs in DNMT3A intron 2:", nrow(cpg_intron2), "\n")

# ============================================================
# Convert to BED coordinates
#
# BED:
#   start = 0-based
#   end   = half-open
#
# For a CpG at 1-based positions x:(x+1):
#   BED = (x-1):(x+1)
# ============================================================

cpg_bed <- cpg |>
  transmute(
    chrom,
    start = start_1based - 1,
    end = end_1based,
    name = paste0("CpG_", row_number())
  )

cpg_intron2_bed <- cpg |>
  filter(in_intron2) |>
  transmute(
    chrom,
    start = start_1based - 1,
    end = end_1based,
    name = paste0("CpG_intron2_", row_number())
  )

# ============================================================
# Save reference annotations
# ============================================================

write_tsv(
  cpg_bed,
  here("metadata", "DNMT3A_window_hg38_CpGs.bed"),
  col_names = FALSE
)

write_tsv(
  cpg_intron2_bed,
  here("metadata", "DNMT3A_intron2_hg38_CpGs.bed"),
  col_names = FALSE
)

cat("\nSaved:\n")
cat("metadata/DNMT3A_window_hg38_CpGs.bed\n")
cat("metadata/DNMT3A_intron2_hg38_CpGs.bed\n")

# ============================================================
# Compare reference CpGs with deposited WGBS intervals
# Representative sample: PBMC 103
# ============================================================

bigwig_manifest <- read_tsv(
  here("metadata", "regional_1d_bigwig_tracks.tsv"),
  show_col_types = FALSE
)

wgbs_103 <- bigwig_manifest |>
  filter(
    track_type == "WGBS",
    group == "PBMC",
    sample_id == "103"
  )

cat("\nPBMC 103 WGBS tracks found:", nrow(wgbs_103), "\n")

# Read the regional BigWig
wgbs_intervals <- cpp11bigwig::read_bigwig(
  wgbs_103$bigwig_file,
  chrom = chrom,
  start = window_start,
  end = window_end
)

cat("WGBS intervals:", nrow(wgbs_intervals), "\n")

# ============================================================
# Count reference CpGs overlapping each WGBS interval
#
# cpp11bigwig intervals use genomic start/end coordinates.
# We compare them using GenomicRanges to avoid assumptions
# about coordinate conventions during overlap calculations.
# ============================================================

library(GenomicRanges)

wgbs_gr <- GRanges(
  seqnames = wgbs_intervals$chrom,
  ranges = IRanges(
    start = wgbs_intervals$start + 1,
    end = wgbs_intervals$end
  )
)

cpg_gr <- GRanges(
  seqnames = cpg$chrom,
  ranges = IRanges(
    start = cpg$start_1based,
    end = cpg$end_1based
  )
)

hits <- findOverlaps(
  wgbs_gr,
  cpg_gr,
  ignore.strand = TRUE
)

cpg_count <- tabulate(
  queryHits(hits),
  nbins = length(wgbs_gr)
)

wgbs_intervals$cpg_count <- cpg_count

cat("\nReference CpGs per deposited WGBS interval:\n")

wgbs_intervals |>
  count(cpg_count) |>
  arrange(cpg_count) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ============================================================
# Restrict comparison to DNMT3A intron 2
# ============================================================

intron2_gr <- GRanges(
  seqnames = chrom,
  ranges = IRanges(
    start = intron2_start,
    end = intron2_end
  )
)

intron2_hits <- findOverlaps(
  wgbs_gr,
  intron2_gr,
  ignore.strand = TRUE
)

wgbs_intron2 <- wgbs_intervals[
  unique(queryHits(intron2_hits)),
]

cat("\nWGBS intervals overlapping intron 2:",
    nrow(wgbs_intron2), "\n")

cat("\nReference CpGs per WGBS interval in intron 2:\n")

wgbs_intron2 |>
  count(cpg_count) |>
  arrange(cpg_count) |>
  as.data.frame() |>
  print(row.names = FALSE)
# ============================================================
# Inspect actual WGBS BigWig interval widths
# ============================================================

wgbs_intervals <- wgbs_intervals |>
  mutate(
    interval_width = end - start
  )

cat("\nWGBS interval-width distribution:\n")

wgbs_intervals |>
  count(interval_width) |>
  arrange(interval_width) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nWGBS interval-width summary:\n")

summary(wgbs_intervals$interval_width)

cat("\nExamples of intervals containing many CpGs:\n")

wgbs_intervals |>
  filter(cpg_count >= 5) |>
  arrange(desc(cpg_count)) |>
  select(chrom, start, end, interval_width, value, cpg_count) |>
  head(20) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ============================================================
# Save WGBS / CpG overlap analysis
# ============================================================

wgbs_cpg_summary <- bind_rows(
  wgbs_intervals |>
    count(cpg_count, name = "n_intervals") |>
    mutate(
      region = "DNMT3A_window",
      percent_intervals = 100 * n_intervals / sum(n_intervals)
    ),
  
  wgbs_intron2 |>
    count(cpg_count, name = "n_intervals") |>
    mutate(
      region = "DNMT3A_intron2",
      percent_intervals = 100 * n_intervals / sum(n_intervals)
    )
) |>
  select(
    region,
    cpg_count,
    n_intervals,
    percent_intervals
  ) |>
  arrange(region, cpg_count)

write_tsv(
  wgbs_cpg_summary,
  here("metadata", "WGBS_PBMC103_CpG_interval_summary.tsv")
)

cat("\nSaved summary table:\n")
cat("metadata/WGBS_PBMC103_CpG_interval_summary.tsv\n")

print(
  wgbs_cpg_summary,
  n = Inf
)

# ============================================================
# Save compact WGBS data-structure summary
# ============================================================

wgbs_structure_summary <- tibble(
  metric = c(
    "Reference CpGs",
    "Stored WGBS BigWig intervals",
    "Stored intervals exactly 10 bp",
    "Stored intervals longer than 10 bp",
    "Percent stored intervals exactly 10 bp",
    "Median stored interval width (bp)",
    "Mean stored interval width (bp)",
    "Maximum stored interval width (bp)"
  ),
  value = c(
    nrow(cpg),
    nrow(wgbs_intervals),
    sum(wgbs_intervals$interval_width == 10),
    sum(wgbs_intervals$interval_width > 10),
    100 * mean(wgbs_intervals$interval_width == 10),
    median(wgbs_intervals$interval_width),
    mean(wgbs_intervals$interval_width),
    max(wgbs_intervals$interval_width)
  )
)

write_tsv(
  wgbs_structure_summary,
  here("metadata", "WGBS_PBMC103_data_structure_summary.tsv")
)

cat("\nWGBS data-structure summary:\n")

wgbs_structure_summary |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nSaved:\n")
cat("metadata/WGBS_PBMC103_data_structure_summary.tsv\n")