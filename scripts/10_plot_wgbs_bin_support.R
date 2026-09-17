library(tidyverse)
library(here)

# ============================================================
# Step 1: load WGBS 10-bp group-summary data
# ============================================================

wgbs <- read_tsv(
  here(
    "metadata",
    "WGBS_group_summary_10bp.tsv"
  ),
  show_col_types = FALSE
)

cat("Rows in complete WGBS summary table:", nrow(wgbs), "\n")

cat("\nColumns:\n")
print(names(wgbs))

cat("\nGroups:\n")
print(table(wgbs$group))

# ============================================================
# Step 2: define genomic regions
# hg38 coordinates, 1-based inclusive
# ============================================================

dnmt3a_start <- 25227874
dnmt3a_end   <- 25342590

intron2_start <- 25300244
intron2_end   <- 25313912

cat("\nRegions:\n")
cat(
  "Whole DNMT3A:",
  paste0("chr2:", dnmt3a_start, "-", dnmt3a_end),
  "\n"
)

cat(
  "DNMT3A intron 2:",
  paste0("chr2:", intron2_start, "-", intron2_end),
  "\n"
)

# ============================================================
# Step 3: select bins overlapping each region
#
# WGBS table uses BED-like coordinates:
# start = 0-based
# end   = end-exclusive
#
# Region coordinates above are 1-based inclusive.
# ============================================================

dnmt3a_start_bed <- dnmt3a_start - 1
intron2_start_bed <- intron2_start - 1

wgbs_dnmt3a <- wgbs |>
  filter(
    chrom == "chr2",
    start < dnmt3a_end,
    end > dnmt3a_start_bed
  )

wgbs_intron2 <- wgbs |>
  filter(
    chrom == "chr2",
    start < intron2_end,
    end > intron2_start_bed
  )

# ============================================================
# Step 4: inspect regional data
# ============================================================

cat("\nWhole DNMT3A — observed group bins:\n")

wgbs_dnmt3a |>
  count(group) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nIntron 2 — observed group bins:\n")

wgbs_intron2 |>
  count(group) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nWhole DNMT3A — sample-support distribution:\n")

wgbs_dnmt3a |>
  count(group, n_samples) |>
  arrange(group, n_samples) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nIntron 2 — sample-support distribution:\n")

wgbs_intron2 |>
  count(group, n_samples) |>
  arrange(group, n_samples) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ============================================================
# Step 5: patient-level WGBS bin coverage
# ============================================================

wgbs_sample <- read_tsv(
  here(
    "metadata",
    "WGBS_sample_level_10bp.tsv"
  ),
  show_col_types = FALSE
)

cat(
  "\nRows in sample-level 10-bp table:",
  nrow(wgbs_sample),
  "\n"
)

cat("\nColumns in sample-level table:\n")
print(names(wgbs_sample))

# ------------------------------------------------------------
# Select bins overlapping the whole DNMT3A genomic span
# ------------------------------------------------------------

sample_dnmt3a <- wgbs_sample |>
  filter(
    chrom == "chr2",
    start < dnmt3a_end,
    end > dnmt3a_start_bed
  )

# ------------------------------------------------------------
# Select bins overlapping DNMT3A intron 2
# ------------------------------------------------------------

sample_intron2 <- wgbs_sample |>
  filter(
    chrom == "chr2",
    start < intron2_end,
    end > intron2_start_bed
  )

# ------------------------------------------------------------
# Count observed 10-bp bins per patient
# ------------------------------------------------------------

patient_coverage_dnmt3a <- sample_dnmt3a |>
  count(
    group,
    sample_id,
    name = "n_observed_bins"
  ) |>
  arrange(group, sample_id)

patient_coverage_intron2 <- sample_intron2 |>
  count(
    group,
    sample_id,
    name = "n_observed_bins"
  ) |>
  arrange(group, sample_id)

cat("\nObserved 10-bp bins per patient — whole DNMT3A:\n")

patient_coverage_dnmt3a |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nObserved 10-bp bins per patient — DNMT3A intron 2:\n")

patient_coverage_intron2 |>
  as.data.frame() |>
  print(row.names = FALSE)

# ============================================================
# Step 6: save patient-level statistics and create plots
# ============================================================

# ------------------------------------------------------------
# Save underlying patient-level statistics
# ------------------------------------------------------------

write_tsv(
  patient_coverage_dnmt3a,
  here(
    "metadata",
    "WGBS_DNMT3A_gene_span_patient_coverage.tsv"
  )
)

write_tsv(
  patient_coverage_intron2,
  here(
    "metadata",
    "WGBS_DNMT3A_intron2_patient_coverage.tsv"
  )
)

cat("\nSaved patient-level statistics:\n")
cat("metadata/WGBS_DNMT3A_gene_span_patient_coverage.tsv\n")
cat("metadata/WGBS_DNMT3A_intron2_patient_coverage.tsv\n")


# ------------------------------------------------------------
# Plotting function
#
# One bar = one individual sample
# Bar height = number of observed WGBS 10-bp bins
# ------------------------------------------------------------

plot_patient_coverage <- function(data, plot_title, plot_subtitle) {
  
  ggplot(
    data,
    aes(
      x = factor(sample_id, levels = unique(sample_id)),
      y = n_observed_bins,
      fill = group
    )
  ) +
    geom_col(
      width = 0.75,
      show.legend = FALSE
    ) +
    geom_text(
      aes(label = n_observed_bins),
      vjust = -0.4,
      size = 3.2
    ) +
    facet_grid(
      . ~ group,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_fill_manual(
      values = c(
        "AML" = "#B22222",
        "PBMC" = "#1E5AB4"
      )
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.10)),
      labels = scales::comma
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Sample",
      y = "Observed WGBS 10-bp bins",
      caption = paste0(
        "GSE152136 processed WGBS; ",
        "bars show observed 10-bp bins per individual sample."
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = element_text(
        size = 11
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 11
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.title = element_text(
        face = "bold"
      ),
      plot.caption = element_text(
        hjust = 0,
        size = 9
      )
    )
}


# ------------------------------------------------------------
# Whole DNMT3A genomic span
# ------------------------------------------------------------

p_dnmt3a_patient_coverage <- plot_patient_coverage(
  patient_coverage_dnmt3a,
  plot_title = "WGBS coverage across the DNMT3A genomic span",
  plot_subtitle = paste0(
    "hg38 chr2:",
    format(dnmt3a_start, big.mark = ","),
    "\u2013",
    format(dnmt3a_end, big.mark = ",")
  )
)

print(p_dnmt3a_patient_coverage)


# ------------------------------------------------------------
# DNMT3A intron 2
# ------------------------------------------------------------

p_intron2_patient_coverage <- plot_patient_coverage(
  patient_coverage_intron2,
  plot_title = "WGBS coverage across DNMT3A intron 2",
  plot_subtitle = paste0(
    "hg38 chr2:",
    format(intron2_start, big.mark = ","),
    "\u2013",
    format(intron2_end, big.mark = ",")
  )
)

print(p_intron2_patient_coverage)


# ------------------------------------------------------------
# Save figures
# ------------------------------------------------------------

ggsave(
  filename = here(
    "figures",
    "WGBS_DNMT3A_gene_span_patient_coverage.png"
  ),
  plot = p_dnmt3a_patient_coverage,
  width = 10,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = here(
    "figures",
    "WGBS_DNMT3A_gene_span_patient_coverage.pdf"
  ),
  plot = p_dnmt3a_patient_coverage,
  width = 10,
  height = 5.5
)

ggsave(
  filename = here(
    "figures",
    "WGBS_DNMT3A_intron2_patient_coverage.png"
  ),
  plot = p_intron2_patient_coverage,
  width = 10,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = here(
    "figures",
    "WGBS_DNMT3A_intron2_patient_coverage.pdf"
  ),
  plot = p_intron2_patient_coverage,
  width = 10,
  height = 5.5
)

cat("\nSaved figures:\n")
cat("figures/WGBS_DNMT3A_gene_span_patient_coverage.png\n")
cat("figures/WGBS_DNMT3A_gene_span_patient_coverage.pdf\n")
cat("figures/WGBS_DNMT3A_intron2_patient_coverage.png\n")
cat("figures/WGBS_DNMT3A_intron2_patient_coverage.pdf\n")

# ============================================================
# Step 7: regional WGBS bin-support statistics
# ============================================================

# ------------------------------------------------------------
# Count bins at each exact sample-support level
# ------------------------------------------------------------

bin_support_dnmt3a <- wgbs_dnmt3a |>
  count(
    group,
    n_samples,
    name = "n_bins"
  ) |>
  mutate(
    passes_threshold = case_when(
      group == "AML"  & n_samples >= 16 ~ TRUE,
      group == "PBMC" & n_samples == 2  ~ TRUE,
      TRUE ~ FALSE
    )
  ) |>
  arrange(group, n_samples)

bin_support_intron2 <- wgbs_intron2 |>
  count(
    group,
    n_samples,
    name = "n_bins"
  ) |>
  mutate(
    passes_threshold = case_when(
      group == "AML"  & n_samples >= 16 ~ TRUE,
      group == "PBMC" & n_samples == 2  ~ TRUE,
      TRUE ~ FALSE
    )
  ) |>
  arrange(group, n_samples)


# ------------------------------------------------------------
# Calculate regional threshold summaries
# ------------------------------------------------------------

support_summary_dnmt3a <- bin_support_dnmt3a |>
  group_by(group) |>
  summarise(
    total_observed_bins = sum(n_bins),
    bins_meeting_threshold = sum(
      n_bins[passes_threshold]
    ),
    percent_meeting_threshold =
      100 * bins_meeting_threshold / total_observed_bins,
    .groups = "drop"
  ) |>
  mutate(
    region = "Whole DNMT3A genomic span",
    threshold = case_when(
      group == "AML" ~ ">=16/18",
      group == "PBMC" ~ "2/2"
    )
  ) |>
  select(
    region,
    group,
    threshold,
    total_observed_bins,
    bins_meeting_threshold,
    percent_meeting_threshold
  )

support_summary_intron2 <- bin_support_intron2 |>
  group_by(group) |>
  summarise(
    total_observed_bins = sum(n_bins),
    bins_meeting_threshold = sum(
      n_bins[passes_threshold]
    ),
    percent_meeting_threshold =
      100 * bins_meeting_threshold / total_observed_bins,
    .groups = "drop"
  ) |>
  mutate(
    region = "DNMT3A intron 2",
    threshold = case_when(
      group == "AML" ~ ">=16/18",
      group == "PBMC" ~ "2/2"
    )
  ) |>
  select(
    region,
    group,
    threshold,
    total_observed_bins,
    bins_meeting_threshold,
    percent_meeting_threshold
  )

support_summary_regions <- bind_rows(
  support_summary_dnmt3a,
  support_summary_intron2
)


# ------------------------------------------------------------
# Save underlying statistics
# ------------------------------------------------------------

write_tsv(
  bin_support_dnmt3a,
  here(
    "metadata",
    "WGBS_DNMT3A_gene_span_bin_support.tsv"
  )
)

write_tsv(
  bin_support_intron2,
  here(
    "metadata",
    "WGBS_DNMT3A_intron2_bin_support.tsv"
  )
)

write_tsv(
  support_summary_regions,
  here(
    "metadata",
    "WGBS_DNMT3A_regional_support_summary.tsv"
  )
)


# ------------------------------------------------------------
# Print threshold summary
# ------------------------------------------------------------

cat("\nRegional WGBS support summary:\n")

support_summary_regions |>
  mutate(
    percent_meeting_threshold =
      round(percent_meeting_threshold, 2)
  ) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\nSaved regional support statistics.\n")

# ============================================================
# Step 8: plot regional WGBS bin-support distributions
# ============================================================

# ------------------------------------------------------------
# Plotting function
#
# X-axis = exact number of samples contributing data to a bin
# Y-axis = number of 10-bp bins with that exact support
#
# Gray = below summary-track threshold
# Red/blue = retained for summary track
# ------------------------------------------------------------

plot_bin_support <- function(data, summary_data,
                             plot_title, plot_subtitle) {
  
  plot_data <- data |>
    mutate(
      support_class = case_when(
        passes_threshold & group == "AML"  ~ "AML retained",
        passes_threshold & group == "PBMC" ~ "PBMC retained",
        TRUE ~ "Below threshold"
      )
    )
  
  annotation_data <- summary_data |>
    mutate(
      label = paste0(
        scales::comma(bins_meeting_threshold),
        "/",
        scales::comma(total_observed_bins),
        " retained\n(",
        sprintf("%.1f", percent_meeting_threshold),
        "%)"
      )
    )
  
  ggplot(
    plot_data,
    aes(
      x = factor(n_samples),
      y = n_bins,
      fill = support_class
    )
  ) +
    geom_col(
      width = 0.75,
      show.legend = FALSE
    ) +
    geom_text(
      aes(label = n_bins),
      vjust = -0.35,
      size = 3.2
    ) +
    facet_grid(
      . ~ group,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_fill_manual(
      values = c(
        "AML retained" = "#B22222",
        "PBMC retained" = "#1E5AB4",
        "Below threshold" = "grey75"
      )
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.12)),
      labels = scales::comma
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Number of samples with data",
      y = "Number of observed WGBS 10-bp bins",
      caption = paste0(
        "Colored bars meet the WGBS summary-track support threshold: ",
        "AML \u226516/18; PBMC 2/2. ",
        "Gray bars fall below the threshold."
      )
    ) +
    geom_text(
      data = annotation_data,
      aes(
        x = Inf,
        y = Inf,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 1.1,
      vjust = 1.2,
      size = 3.2,
      fontface = "bold",
      lineheight = 0.9
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = element_text(
        size = 11
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 11
      ),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5
      ),
      axis.title = element_text(
        face = "bold"
      ),
      plot.caption = element_text(
        hjust = 0,
        size = 9
      )
    )
}


# ------------------------------------------------------------
# Whole DNMT3A genomic span
# ------------------------------------------------------------

p_dnmt3a_bin_support <- plot_bin_support(
  bin_support_dnmt3a,
  support_summary_dnmt3a,
  plot_title = "WGBS bin support across the DNMT3A genomic span",
  plot_subtitle = paste0(
    "hg38 chr2:",
    format(dnmt3a_start, big.mark = ","),
    "\u2013",
    format(dnmt3a_end, big.mark = ",")
  )
)

print(p_dnmt3a_bin_support)


# ------------------------------------------------------------
# DNMT3A intron 2
# ------------------------------------------------------------

p_intron2_bin_support <- plot_bin_support(
  bin_support_intron2,
  support_summary_intron2,
  plot_title = "WGBS bin support across DNMT3A intron 2",
  plot_subtitle = paste0(
    "hg38 chr2:",
    format(intron2_start, big.mark = ","),
    "\u2013",
    format(intron2_end, big.mark = ",")
  )
)

print(p_intron2_bin_support)


# ------------------------------------------------------------
# Save figures
# ------------------------------------------------------------

ggsave(
  here(
    "figures",
    "WGBS_DNMT3A_gene_span_bin_support.png"
  ),
  p_dnmt3a_bin_support,
  width = 10,
  height = 5.5,
  dpi = 300
)

ggsave(
  here(
    "figures",
    "WGBS_DNMT3A_gene_span_bin_support.pdf"
  ),
  p_dnmt3a_bin_support,
  width = 10,
  height = 5.5
)

ggsave(
  here(
    "figures",
    "WGBS_DNMT3A_intron2_bin_support.png"
  ),
  p_intron2_bin_support,
  width = 10,
  height = 5.5,
  dpi = 300
)

ggsave(
  here(
    "figures",
    "WGBS_DNMT3A_intron2_bin_support.pdf"
  ),
  p_intron2_bin_support,
  width = 10,
  height = 5.5
)

cat("\nSaved bin-support figures:\n")
cat("figures/WGBS_DNMT3A_gene_span_bin_support.png\n")
cat("figures/WGBS_DNMT3A_gene_span_bin_support.pdf\n")
cat("figures/WGBS_DNMT3A_intron2_bin_support.png\n")
cat("figures/WGBS_DNMT3A_intron2_bin_support.pdf\n")