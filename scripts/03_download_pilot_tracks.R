# ============================================================
# GSE152136 DNMT3A multi-omics project
# Step 3: Pilot download of browser-ready 1D tracks
# ============================================================

library(tidyverse)
library(here)
library(rtracklayer)

# Load processed 1D track manifest
tracks_1d <- read_tsv(
  here("metadata", "processed_tracks_1d.tsv"),
  show_col_types = FALSE
)

# Inspect dimensions without tibble-printing issues
cat(
  "Number of 1D tracks:", nrow(tracks_1d),
  "\n"
)

# ============================================================
# Select pilot tracks: PBMC sample 103
# ============================================================

pilot_tracks <- tracks_1d |>
  filter(sample_id == "103") |>
  filter(
    assay == "WGBS" |
      assay == "ATAC-seq" |
      (assay == "CUT&Tag" &
         mark %in% c("H3K27ac", "H3K27me3", "CTCF"))
  )

cat(
  "Pilot tracks selected:", nrow(pilot_tracks),
  "\n"
)

as.data.frame(
  pilot_tracks |>
    select(GSM, sample_id, group, assay, mark, fname)
)

# ============================================================
# Check remote file sizes BEFORE downloading
# ============================================================

get_remote_size <- function(url) {
  
  h <- httr2::request(url) |>
    httr2::req_method("HEAD") |>
    httr2::req_perform()
  
  size_bytes <- httr2::resp_header(h, "content-length")
  
  if (is.null(size_bytes)) {
    return(NA_real_)
  }
  
  as.numeric(size_bytes)
}

pilot_sizes <- pilot_tracks |>
  mutate(
    size_bytes = purrr::map_dbl(url, get_remote_size),
    size_MB = size_bytes / 1024^2,
    size_GB = size_bytes / 1024^3
  )

as.data.frame(
  pilot_sizes |>
    select(
      sample_id,
      assay,
      mark,
      fname,
      size_MB,
      size_GB
    )
)
cat(
  "\nTotal pilot download:",
  round(sum(pilot_sizes$size_GB, na.rm = TRUE), 3),
  "GB\n"
)

# ============================================================
# Download pilot bigWig files
# Measure download time for each file
# ============================================================

data_root <- "C:/Genomics/GSE152136"

pilot_dir <- file.path(
  data_root,
  "pilot",
  "PBMC_103"
)

dir.create(
  pilot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "\nPilot download directory:",
  pilot_dir,
  "\n"
)

# Start total timer
total_start <- Sys.time()

download_times <- numeric(nrow(pilot_tracks))

for (i in seq_len(nrow(pilot_tracks))) {
  
  file_url <- pilot_tracks$url[i]
  file_name <- pilot_tracks$fname[i]
  
  destination <- file.path(
    pilot_dir,
    file_name
  )
  
  if (file.exists(destination)) {
    
    message(
      "Already exists, skipping: ",
      file_name
    )
    
    download_times[i] <- NA_real_
    
  } else {
    
    message(
      "Downloading: ",
      file_name
    )
    
    file_start <- Sys.time()
    
    download.file(
      url = file_url,
      destfile = destination,
      mode = "wb",
      method = "libcurl"
    )
    
    file_end <- Sys.time()
    
    download_times[i] <- as.numeric(
      difftime(
        file_end,
        file_start,
        units = "secs"
      )
    )
    
    message(
      "Finished in ",
      round(download_times[i], 1),
      " seconds."
    )
  }
}

# Stop total timer
total_end <- Sys.time()

total_seconds <- as.numeric(
  difftime(
    total_end,
    total_start,
    units = "secs"
  )
)

cat(
  "\nPilot download finished.",
  "\nTotal elapsed time:",
  round(total_seconds, 1),
  "seconds",
  "\nTotal elapsed time:",
  round(total_seconds / 60, 2),
  "minutes\n"
)

# ============================================================
# Verify downloaded pilot files
# ============================================================

pilot_downloads <- tibble(
  file = list.files(
    pilot_dir,
    full.names = TRUE
  )
) |>
  mutate(
    filename = basename(file),
    size_bytes = file.info(file)$size,
    size_MB = size_bytes / 1024^2
  )

cat(
  "\nDownloaded files:", nrow(pilot_downloads),
  "\nTotal size:",
  round(sum(pilot_downloads$size_MB), 1),
  "MB\n"
)

as.data.frame(
  pilot_downloads |>
    select(filename, size_MB)
)

# ============================================================
# Test bigWig access at the DNMT3A locus
# hg38
# ============================================================

library(GenomicRanges)
library(rtracklayer)

dnmt3a_region <- GRanges(
  seqnames = "chr2",
  ranges = IRanges(
    start = 25227874,
    end   = 25341925
  )
)

dnmt3a_region

dnmt3a_imports <- lapply(
  pilot_downloads$file,
  function(bigwig_file) {
    
    message(
      "Reading DNMT3A from: ",
      basename(bigwig_file)
    )
    
    rtracklayer::import(
      bigwig_file,
      format = "BigWig",
      which = dnmt3a_region
    )
  }
)

names(dnmt3a_imports) <- pilot_downloads$filename

# ============================================================
# DNMT3A regional window
# hg38: DNMT3A locus plus ~500 kb flanking sequence
# ============================================================

dnmt3a_window_bed <- here(
  "metadata",
  "DNMT3A_window_hg38.bed"
)

writeLines(
  "chr2\t24727873\t25842590\tDNMT3A_window",
  dnmt3a_window_bed
)

cat(
  "DNMT3A window BED:\n",
  readLines(dnmt3a_window_bed),
  "\n"
)

remote_test_url <- pilot_tracks |>
  filter(
    sample_id == "103",
    assay == "WGBS"
  ) |>
  pull(url)

cat(
  "Remote test URL obtained:",
  length(remote_test_url) == 1,
  "\n"
)

cat("BED used for test:\n")
cat(readLines(dnmt3a_window_bed), sep = "\n")
cat("\n\n")

remote_start <- Sys.time()

remote_test <- megadepth::get_coverage(
  remote_test_url,
  annotation = dnmt3a_window_bed,
  op = "mean"
)

remote_end <- Sys.time()

cat(
  "\nElapsed:",
  round(
    as.numeric(difftime(remote_end, remote_start, units = "secs")),
    2
  ),
  "seconds\n"
)

cat(
  "remote_test exists:",
  exists("remote_test"),
  "\n"
)
# ============================================================
# Extract regional WGBS signal from remote GEO BigWig
# PBMC 103, DNMT3A +/- 500 kb
# ============================================================

regional_dir <- "C:/Genomics/GSE152136/regional/PBMC_103"

dir.create(
  regional_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

regional_prefix <- file.path(
  regional_dir,
  "PBMC_103_WGBS_DNMT3A"
)

extract_start <- Sys.time()

megadepth::megadepth_cmd(
  paste0(
    '"', remote_test_url, '" ',
    '--annotation "', dnmt3a_window_bed, '" ',
    '--prefix "', regional_prefix, '"'
  )
)

extract_end <- Sys.time()

cat(
  "\nExtraction time:",
  round(
    as.numeric(difftime(extract_end, extract_start, units = "secs")),
    2
  ),
  "seconds\n"
)

cat("\nFiles created:\n")

print(
  list.files(
    regional_dir,
    full.names = TRUE
  )
)

# ============================================================
# Small WGBS resolution test around DNMT3A intron 2
# NM_022552.5, hg38
# ============================================================

test_10bp_bed <- here(
  "metadata",
  "DNMT3A_intron2_test_10bp.bed"
)

test_starts <- seq(
  from = 25300243,
  to   = 25300432,
  by   = 10
)

test_windows <- tibble(
  chr   = "chr2",
  start = test_starts,
  end   = test_starts + 10
)

write.table(
  test_windows,
  test_10bp_bed,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat(readLines(test_10bp_bed), sep = "\n")

test_presence_bed <- here(
  "metadata",
  "WGBS_presence_test.bed"
)

writeLines(
  c(
    "chr2\t25300243\t25300253\tpositive_test",
    "chr2\t25300253\t25300263\tzero_test"
  ),
  test_presence_bed
)

cat(readLines(test_presence_bed), sep = "\n")

cat("\n--- MIN ---\n")

megadepth::megadepth_cmd(
  paste0(
    '"', remote_test_url, '" ',
    '--annotation "', test_presence_bed, '" ',
    '--op min'
  )
)

cat("\n--- MAX ---\n")

megadepth::megadepth_cmd(
  paste0(
    '"', remote_test_url, '" ',
    '--annotation "', test_presence_bed, '" ',
    '--op max'
  )
)

# ============================================================
# Test cpp11bigwig directly on remote GEO WGBS BigWig
# PBMC 103, DNMT3A +/- 500 kb
# ============================================================

remote_cpp_start <- Sys.time()

wgbs_region_remote <- cpp11bigwig::read_bigwig(
  remote_test_url,
  chrom = "chr2",
  start = 24727873,
  end = 25842590
)

remote_cpp_end <- Sys.time()

cat(
  "\nRemote cpp11bigwig query time:",
  round(
    as.numeric(
      difftime(
        remote_cpp_end,
        remote_cpp_start,
        units = "secs"
      )
    ),
    2
  ),
  "seconds\n"
)

cat(
  "Number of returned intervals:",
  nrow(wgbs_region_remote),
  "\n"
)

# ============================================================
# Remote regional extraction test:
# ATAC-seq + CUT&Tag, PBMC 103
# DNMT3A +/- 500 kb
# ============================================================

test_signal_tracks <- pilot_tracks |>
  filter(
    sample_id == "103",
    assay == "ATAC-seq" |
      (assay == "CUT&Tag" &
         mark %in% c("H3K27ac", "H3K27me3", "CTCF"))
  ) |>
  select(assay, mark, url)

print(as.data.frame(test_signal_tracks))

signal_results <- vector(
  "list",
  nrow(test_signal_tracks)
)

for (i in seq_len(nrow(test_signal_tracks))) {
  
  cat(
    "\nReading",
    test_signal_tracks$assay[i],
    test_signal_tracks$mark[i],
    "...\n"
  )
  
  t0 <- Sys.time()
  
  x <- cpp11bigwig::read_bigwig(
    test_signal_tracks$url[i],
    chrom = "chr2",
    start = 24727873,
    end = 25842590
  )
  
  elapsed <- as.numeric(
    difftime(Sys.time(), t0, units = "secs")
  )
  
  signal_results[[i]] <- x
  
  cat(
    "Intervals:",
    nrow(x),
    "| Time:",
    round(elapsed, 2),
    "seconds\n"
  )
}

# ============================================================
# Hi-C inventory: inspect 32 remote .mcool files
# No files are downloaded in this step
# ============================================================

hic_tracks <- read_tsv(
  here("metadata", "processed_tracks_hic.tsv"),
  show_col_types = FALSE
)

cat(
  "Number of Hi-C tracks:",
  nrow(hic_tracks),
  "\n"
)

print(
  as.data.frame(
    hic_tracks |>
      count(group, name = "n")
  )
)

# ============================================================
# Get remote .mcool file sizes
# ============================================================

hic_sizes <- hic_tracks |>
  mutate(
    size_bytes = purrr::map_dbl(url, get_remote_size),
    size_MB = size_bytes / 1024^2,
    size_GB = size_bytes / 1024^3
  )

cat(
  "\nTotal size of all 32 complete mcool files:",
  round(sum(hic_sizes$size_GB, na.rm = TRUE), 2),
  "GB\n\n"
)

print(
  as.data.frame(
    hic_sizes |>
      group_by(group) |>
      summarise(
        n = n(),
        total_GB = sum(size_GB, na.rm = TRUE),
        mean_GB = mean(size_GB, na.rm = TRUE),
        min_GB = min(size_GB, na.rm = TRUE),
        max_GB = max(size_GB, na.rm = TRUE),
        .groups = "drop"
      )
  )
)

# ============================================================
# Select PBMC 103 Hi-C pilot file
# ============================================================

hic_103 <- hic_sizes |>
  filter(
    sample_id == "103",
    group == "PBMC"
  )

cat(
  "PBMC 103 Hi-C URL:\n",
  hic_103$url,
  "\n\n"
)

cat(
  "Full file size:",
  round(hic_103$size_GB, 3),
  "GB\n"
)

# ============================================================
# Download one Hi-C pilot: PBMC 103
# ============================================================

hic_pilot_dir <- file.path(
  "C:/Genomics/GSE152136/pilot",
  "PBMC_103_HiC"
)

dir.create(
  hic_pilot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

hic_103_file <- file.path(
  hic_pilot_dir,
  hic_103$fname
)

cat(
  "Downloading to:\n",
  hic_103_file,
  "\n\n"
)

t0 <- Sys.time()

download.file(
  hic_103$url,
  destfile = hic_103_file,
  mode = "wb"
)

elapsed <- as.numeric(
  difftime(Sys.time(), t0, units = "mins")
)

cat(
  "\nDownload finished.\n",
  "Elapsed time:",
  round(elapsed, 2),
  "minutes\n"
)

cat(
  "Downloaded size:",
  round(file.info(hic_103_file)$size / 1024^3, 3),
  "GB\n"
)

# ============================================================
# Inspect PBMC 103 .mcool structure
# ============================================================

library(rhdf5)

hic_structure <- h5ls(
  hic_103_file,
  recursive = TRUE
)

cat(
  "Number of HDF5 objects:",
  nrow(hic_structure),
  "\n\n"
)

print(
  as.data.frame(
    hic_structure |>
      filter(
        grepl("resolutions", group) |
          grepl("resolutions", name)
      ) |>
      select(group, name, otype)
  )
)
resolution_groups <- hic_structure |>
  filter(
    group == "/resolutions",
    otype == "H5I_GROUP"
  ) |>
  pull(name)

cat(
  "Available Hi-C resolutions:\n",
  paste(resolution_groups, collapse = "\n"),
  "\n"
)

# ============================================================
# Inspect the 10-kb Hi-C resolution
# ============================================================

hic_10k_structure <- hic_structure |>
  filter(
    grepl("^/resolutions/10000", group) |
      (group == "/resolutions" & name == "10000")
  )

print(
  as.data.frame(
    hic_10k_structure |>
      select(group, name, otype, dclass, dim)
  )
)
# ============================================================
# Inspect attributes of the 10-kb Cooler matrix
# ============================================================

hic_10k_attributes <- h5readAttributes(
  hic_103_file,
  "/resolutions/10000"
)

print(hic_10k_attributes)


# ============================================================
# Locate DNMT3A window in the 10-kb Hi-C matrix
# ============================================================

hic_103_file <- file.path(
  "C:/Genomics/GSE152136/pilot",
  "PBMC_103_HiC",
  "GSM4604272_103.iced.mcool"
)

cat(
  "Hi-C pilot file exists:",
  file.exists(hic_103_file),
  "\n"
)


hic_resolution <- 10000
hic_group <- "/resolutions/10000"

region_chrom <- "chr2"
region_start <- 24727873
region_end   <- 25842590

chrom_names <- h5read(
  hic_103_file,
  paste0(hic_group, "/chroms/name")
)

chrom_offsets <- h5read(
  hic_103_file,
  paste0(hic_group, "/indexes/chrom_offset")
)

chr2_index <- match(region_chrom, chrom_names)

cat(
  "chr2 chromosome index:",
  chr2_index,
  "\n"
)

cat(
  "chr2 global bin range:",
  chrom_offsets[chr2_index],
  "to",
  chrom_offsets[chr2_index + 1] - 1,
  "\n"
)

# ============================================================
# Read only chr2 bin coordinates
# ============================================================

chr2_first_bin <- chrom_offsets[chr2_index]
chr2_last_bin  <- chrom_offsets[chr2_index + 1] - 1

# Cooler bin IDs are 0-based.
# R indices are 1-based, hence the +1.
chr2_r_index <- (chr2_first_bin:chr2_last_bin) + 1

chr2_starts <- h5read(
  hic_103_file,
  paste0(hic_group, "/bins/start"),
  index = list(chr2_r_index)
)

chr2_ends <- h5read(
  hic_103_file,
  paste0(hic_group, "/bins/end"),
  index = list(chr2_r_index)
)

chr2_bins <- tibble(
  bin_id = chr2_first_bin:chr2_last_bin,
  chrom = "chr2",
  start = chr2_starts,
  end = chr2_ends
)

dnmt3a_bins <- chr2_bins |>
  filter(
    end > region_start,
    start < region_end
  )

cat(
  "\nNumber of 10-kb bins overlapping DNMT3A window:",
  nrow(dnmt3a_bins),
  "\n"
)

cat(
  "First global bin ID:",
  min(dnmt3a_bins$bin_id),
  "\n"
)

cat(
  "Last global bin ID:",
  max(dnmt3a_bins$bin_id),
  "\n\n"
)

print(
  as.data.frame(
    bind_rows(
      head(dnmt3a_bins, 3),
      tail(dnmt3a_bins, 3)
    )
  )
)

# ============================================================
# Extract Hi-C pixels for the DNMT3A 10-kb region
# ============================================================

region_first_bin <- min(dnmt3a_bins$bin_id)
region_last_bin  <- max(dnmt3a_bins$bin_id)

# Cooler bin IDs and pixel offsets are 0-based.
# R/HDF5 indices are 1-based.
pixel_offsets <- h5read(
  hic_103_file,
  paste0(hic_group, "/indexes/bin1_offset"),
  index = list(
    c(
      region_first_bin + 1,
      region_last_bin + 2
    )
  )
)

pixel_start <- pixel_offsets[1]
pixel_end   <- pixel_offsets[2]

cat(
  "Pixel offset range:",
  pixel_start,
  "to",
  pixel_end,
  "\n"
)

cat(
  "Pixels to read before regional bin2 filtering:",
  pixel_end - pixel_start,
  "\n"
)

# ============================================================
# Read only those pixels from the HDF5 file
# ============================================================

pixel_r_index <- (pixel_start + 1):pixel_end

region_bin1 <- h5read(
  hic_103_file,
  paste0(hic_group, "/pixels/bin1_id"),
  index = list(pixel_r_index)
)

region_bin2 <- h5read(
  hic_103_file,
  paste0(hic_group, "/pixels/bin2_id"),
  index = list(pixel_r_index)
)

region_count <- h5read(
  hic_103_file,
  paste0(hic_group, "/pixels/count"),
  index = list(pixel_r_index)
)

region_pixels <- tibble(
  bin1_id = region_bin1,
  bin2_id = region_bin2,
  count = region_count
) |>
  filter(
    bin2_id >= region_first_bin,
    bin2_id <= region_last_bin
  )

cat(
  "\nRegional DNMT3A Hi-C pixels:",
  nrow(region_pixels),
  "\n"
)

print(
  head(
    as.data.frame(region_pixels),
    10
  )
)

# ============================================================
# Read balancing weights for the 113 DNMT3A bins
# ============================================================

region_bin_ids <- dnmt3a_bins$bin_id

region_weights <- h5read(
  hic_103_file,
  paste0(hic_group, "/bins/weight"),
  index = list(region_bin_ids + 1)
)

dnmt3a_bins <- dnmt3a_bins |>
  mutate(
    weight = region_weights
  )

cat(
  "Number of bins:",
  nrow(dnmt3a_bins),
  "\n"
)

cat(
  "Bins with missing weights:",
  sum(is.na(dnmt3a_bins$weight)),
  "\n"
)

print(
  head(
    as.data.frame(dnmt3a_bins),
    10
  )
)
# ============================================================
# Add ICE-balanced contact values to regional pixels
# ============================================================

weight_lookup <- setNames(
  dnmt3a_bins$weight,
  dnmt3a_bins$bin_id
)

region_pixels_balanced <- region_pixels |>
  mutate(
    weight1 = weight_lookup[as.character(bin1_id)],
    weight2 = weight_lookup[as.character(bin2_id)],
    balanced = count * weight1 * weight2
  )

cat(
  "Regional pixels:",
  nrow(region_pixels_balanced),
  "\n"
)

cat(
  "Pixels with missing balanced value:",
  sum(is.na(region_pixels_balanced$balanced)),
  "\n\n"
)

print(
  head(
    as.data.frame(region_pixels_balanced),
    10
  )
)

# ============================================================
# Build full symmetric 113 x 113 balanced Hi-C matrix
# ============================================================

n_bins <- nrow(dnmt3a_bins)

bin_index <- setNames(
  seq_len(n_bins),
  dnmt3a_bins$bin_id
)

hic_matrix_balanced <- matrix(
  0,
  nrow = n_bins,
  ncol = n_bins
)

for (i in seq_len(nrow(region_pixels_balanced))) {
  
  r <- bin_index[as.character(region_pixels_balanced$bin1_id[i])]
  c <- bin_index[as.character(region_pixels_balanced$bin2_id[i])]
  
  value <- region_pixels_balanced$balanced[i]
  
  hic_matrix_balanced[r, c] <- value
  hic_matrix_balanced[c, r] <- value
}

cat(
  "Matrix dimensions:",
  nrow(hic_matrix_balanced),
  "x",
  ncol(hic_matrix_balanced),
  "\n"
)

cat(
  "Non-zero matrix cells:",
  sum(hic_matrix_balanced > 0),
  "\n"
)
# ============================================================
# First diagnostic Hi-C heatmap
# ============================================================

image(
  log1p(hic_matrix_balanced),
  axes = FALSE,
  main = "PBMC 103 Hi-C — DNMT3A region — 10 kb"
)


# ============================================================
# Add genomic coordinates to the 113 regional Hi-C bins
# ============================================================

bin_centers <- (
  dnmt3a_bins$start +
    dnmt3a_bins$end
) / 2

cat(
  "First bin center:",
  bin_centers[1],
  "\n"
)

cat(
  "Last bin center:",
  tail(bin_centers, 1),
  "\n"
)

# ============================================================
# DNMT3A annotation coordinates, hg38
# ============================================================

dnmt3a_gene_start <- 25227874
dnmt3a_gene_end   <- 25342590

# NM_022552.5 / ENST00000321117.10
# exact intron 2, 1-based inclusive
intron2_start <- 25300244
intron2_end   <- 25313912

# ============================================================
# Coordinate-labelled balanced Hi-C heatmap
# ============================================================

plot_matrix <- log1p(hic_matrix_balanced)

image(
  x = bin_centers,
  y = bin_centers,
  z = plot_matrix,
  xlab = "chr2 genomic position (bp)",
  ylab = "chr2 genomic position (bp)",
  main = "PBMC 103 Hi-C — DNMT3A region — 10 kb",
  useRaster = TRUE
)

# DNMT3A gene boundaries
abline(
  v = c(dnmt3a_gene_start, dnmt3a_gene_end),
  lty = 2,
  lwd = 1.5
)

abline(
  h = c(dnmt3a_gene_start, dnmt3a_gene_end),
  lty = 2,
  lwd = 1.5
)

# Intron 2 boundaries
abline(
  v = c(intron2_start, intron2_end),
  lty = 3,
  lwd = 2
)

abline(
  h = c(intron2_start, intron2_end),
  lty = 3,
  lwd = 2
)

# ============================================================
# Which 10-kb Hi-C bins overlap DNMT3A and intron 2?
# ============================================================

gene_bins <- dnmt3a_bins |>
  filter(
    end > dnmt3a_gene_start,
    start < dnmt3a_gene_end
  )

intron2_bins <- dnmt3a_bins |>
  filter(
    end > intron2_start,
    start < intron2_end
  )

cat(
  "10-kb bins overlapping DNMT3A gene:",
  nrow(gene_bins),
  "\n"
)

cat(
  "10-kb bins overlapping intron 2:",
  nrow(intron2_bins),
  "\n\n"
)

print(
  as.data.frame(intron2_bins)
)
# ============================================================
# Save PBMC 103 regional Hi-C pilot
# ============================================================

hic_regional_dir <- file.path(
  "C:/Genomics/GSE152136/regional",
  "PBMC_103",
  "HiC"
)

dir.create(
  hic_regional_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Convert HDF5-derived columns to ordinary R vectors
dnmt3a_bins_save <- dnmt3a_bins |>
  transmute(
    bin_id = as.integer(bin_id),
    chrom = as.character(chrom),
    start = as.integer(start),
    end = as.integer(end),
    weight = as.numeric(weight)
  )

region_pixels_save <- region_pixels_balanced |>
  transmute(
    bin1_id = as.integer(bin1_id),
    bin2_id = as.integer(bin2_id),
    count = as.integer(count),
    weight1 = as.numeric(weight1),
    weight2 = as.numeric(weight2),
    balanced = as.numeric(balanced)
  )

# Save the 113 genomic bins
write_tsv(
  dnmt3a_bins_save,
  file.path(
    hic_regional_dir,
    "PBMC_103_HiC_10kb_DNMT3A_bins.tsv"
  )
)

# Save sparse upper-triangle contacts
write_tsv(
  region_pixels_save,
  file.path(
    hic_regional_dir,
    "PBMC_103_HiC_10kb_DNMT3A_pixels.tsv"
  )
)

# Save full symmetric balanced matrix
saveRDS(
  hic_matrix_balanced,
  file.path(
    hic_regional_dir,
    "PBMC_103_HiC_10kb_DNMT3A_balanced_matrix.rds"
  )
)

cat(
  "Regional Hi-C pilot saved to:\n",
  hic_regional_dir,
  "\n\n"
)

print(
  file.info(
    list.files(
      hic_regional_dir,
      full.names = TRUE
    )
  )[, "size", drop = FALSE]
)