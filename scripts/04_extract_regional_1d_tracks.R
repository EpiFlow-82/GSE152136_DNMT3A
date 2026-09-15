library(tidyverse)
library(here)
library(cpp11bigwig)

tracks_1d <- read_tsv(
  here("metadata", "processed_tracks_1d.tsv"),
  show_col_types = FALSE
)

cat(
  "Number of 1D tracks:",
  nrow(tracks_1d),
  "\n\n"
)

cat(
  "Column names:\n"
)

print(
  names(tracks_1d)
)

cat(
  "\nFirst 6 rows:\n"
)

print(
  as.data.frame(
    head(tracks_1d)
  )
)

# ============================================================
# DNMT3A regional coordinates
# hg38, BED coordinates: 0-based, half-open
# ============================================================

region_chrom <- "chr2"
region_start <- 24727873
region_end   <- 25842590

regional_root <- here("tracks", "regional_1d")

dir.create(
  regional_root,
  recursive = TRUE,
  showWarnings = FALSE
)

# Build labels and output paths
tracks_extract <- tracks_1d |>
  mutate(
    track_type = case_when(
      assay == "WGBS" ~ "WGBS",
      assay == "ATAC" ~ "ATAC",
      assay == "CUT&Tag" & !is.na(mark) ~ mark,
      TRUE ~ assay
    ),
    output_dir = file.path(
      regional_root,
      track_type
    ),
    output_file = file.path(
      output_dir,
      paste0(
        GSM, "_",
        sample_id, "_",
        track_type,
        "_DNMT3A.bedGraph"
      )
    )
  )

cat(
  "Tracks prepared for extraction:",
  nrow(tracks_extract),
  "\n\n"
)

print(
  as.data.frame(
    tracks_extract |>
      count(track_type, group)
  )
)


# ============================================================
# Function to extract one regional BigWig track
# ============================================================

extract_one_track <- function(i) {
  
  x <- tracks_extract[i, ]
  
  dir.create(
    x$output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  cat(
    "\n--------------------------------------------------\n",
    "[", i, "/", nrow(tracks_extract), "] ",
    x$GSM, " | ",
    x$sample_id, " | ",
    x$group, " | ",
    x$track_type,
    "\n",
    sep = ""
  )
  
  # Skip files already created successfully
  if (
    file.exists(x$output_file) &&
    file.info(x$output_file)$size > 0
  ) {
    
    cat("Already exists — skipping\n")
    
    return(
      tibble(
        GSM = x$GSM,
        sample_id = x$sample_id,
        group = x$group,
        assay = x$assay,
        mark = x$mark,
        track_type = x$track_type,
        status = "already_exists",
        n_intervals = NA_integer_,
        file_size_bytes = file.info(x$output_file)$size,
        elapsed_seconds = NA_real_,
        output_file = x$output_file,
        error_message = NA_character_
      )
    )
  }
  
  start_time <- Sys.time()
  
  result <- tryCatch({
    
    bw <- cpp11bigwig::read_bigwig(
      x$url,
      chrom = region_chrom,
      start = region_start,
      end = region_end
    )
    
    # Convert any array/matrix columns returned by HDF5/C++
    # infrastructure into ordinary R vectors.
    bw <- tibble(
      chrom = as.character(bw$chrom),
      start = as.integer(bw$start),
      end = as.integer(bw$end),
      value = as.numeric(bw$value)
    )
    
    write_tsv(
      bw,
      x$output_file,
      col_names = FALSE
    )
    
    elapsed <- as.numeric(
      difftime(
        Sys.time(),
        start_time,
        units = "secs"
      )
    )
    
    cat(
      "Intervals:", nrow(bw),
      "| Size:",
      round(file.info(x$output_file)$size / 1024^2, 3),
      "MB",
      "| Time:",
      round(elapsed, 2),
      "sec\n"
    )
    
    tibble(
      GSM = x$GSM,
      sample_id = x$sample_id,
      group = x$group,
      assay = x$assay,
      mark = x$mark,
      track_type = x$track_type,
      status = "success",
      n_intervals = nrow(bw),
      file_size_bytes = file.info(x$output_file)$size,
      elapsed_seconds = elapsed,
      output_file = x$output_file,
      error_message = NA_character_
    )
    
  }, error = function(e) {
    
    elapsed <- as.numeric(
      difftime(
        Sys.time(),
        start_time,
        units = "secs"
      )
    )
    
    cat(
      "ERROR:",
      conditionMessage(e),
      "\n"
    )
    
    tibble(
      GSM = x$GSM,
      sample_id = x$sample_id,
      group = x$group,
      assay = x$assay,
      mark = x$mark,
      track_type = x$track_type,
      status = "error",
      n_intervals = NA_integer_,
      file_size_bytes = NA_real_,
      elapsed_seconds = elapsed,
      output_file = x$output_file,
      error_message = conditionMessage(e)
    )
  })
  
  result
}




# ============================================================
# Three-track extraction test
# ============================================================

test_indices <- c(
  which(tracks_extract$track_type == "WGBS")[1],
  which(tracks_extract$track_type == "ATAC-seq")[1],
  which(tracks_extract$track_type == "H3K27ac")[1]
)

test_results <- map_dfr(
  test_indices,
  extract_one_track
)

print(
  as.data.frame(test_results)
)

# ============================================================
# Extract all 69 regional 1D tracks
# ============================================================

all_results <- map_dfr(
  seq_len(nrow(tracks_extract)),
  extract_one_track
)

cat(
  "\nExtraction finished.\n\n"
)

print(
  as.data.frame(
    all_results |>
      count(status)
  )
)

cat(
  "\nTotal regional track size:\n"
)

total_bytes <- sum(
  all_results$file_size_bytes,
  na.rm = TRUE
)

cat(
  round(total_bytes / 1024^2, 2),
  "MB\n"
)

# ============================================================
# Save regional 1D extraction manifest
# ============================================================

regional_manifest <- all_results |>
  left_join(
    tracks_extract |>
      select(
        GSM,
        fname,
        url,
        GSE,
        title
      ),
    by = "GSM"
  ) |>
  arrange(
    track_type,
    group,
    sample_id
  )

write_tsv(
  regional_manifest,
  here(
    "metadata",
    "regional_1d_tracks.tsv"
  )
)

cat(
  "Regional 1D manifest saved:\n",
  here("metadata", "regional_1d_tracks.tsv"),
  "\n\n"
)

cat(
  "Tracks in manifest:",
  nrow(regional_manifest),
  "\n"
)

cat(
  "Tracks with extraction errors:",
  sum(regional_manifest$status == "error"),
  "\n"
)