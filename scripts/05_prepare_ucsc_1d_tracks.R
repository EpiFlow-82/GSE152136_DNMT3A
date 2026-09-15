print(
  Sys.which("bedGraphToBigWig")
)
print(
  Sys.which("python")
)

print(
  Sys.which("py")
)

cat("Python version:\n")

print(
  system2(
    Sys.which("python"),
    "--version",
    stdout = TRUE,
    stderr = TRUE
  )
)

cat("\nPip version:\n")

print(
  system2(
    Sys.which("python"),
    c("-m", "pip", "--version"),
    stdout = TRUE,
    stderr = TRUE
  )
)

cat("WSL executable:\n")
print(
  Sys.which("wsl")
)

cat("\nWSL status:\n")
print(
  system2(
    "wsl",
    "--status",
    stdout = TRUE,
    stderr = TRUE
  )
)
cat("\nInstalled WSL distributions:\n")

print(
  system2(
    "wsl",
    c("--list", "--verbose"),
    stdout = TRUE,
    stderr = TRUE
  )
)
cat("WSL distributions after restart:\n")

print(
  system2(
    "wsl",
    c("--list", "--verbose"),
    stdout = TRUE,
    stderr = TRUE
  )
)

library(tidyverse)
library(here)

regional_manifest <- read_tsv(
  here("metadata", "regional_1d_tracks.tsv"),
  show_col_types = FALSE
)

windows_to_wsl <- function(path) {
  path <- normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )
  
  drive <- substr(path, 1, 1)
  rest <- substr(path, 3, nchar(path))
  
  paste0(
    "/mnt/",
    tolower(drive),
    rest
  )
}

convert_one_bigwig <- function(i) {
  
  x <- regional_manifest[i, ]
  
  bedgraph_win <- x$output_file
  
  bigwig_win <- sub(
    "\\.bedGraph$",
    ".bw",
    bedgraph_win
  )
  
  bedgraph_wsl <- windows_to_wsl(bedgraph_win)
  bigwig_wsl <- windows_to_wsl(bigwig_win)
  
  cat(
    "\n[", i, "/", nrow(regional_manifest), "] ",
    x$GSM, " | ",
    x$sample_id, " | ",
    x$track_type,
    "\n",
    sep = ""
  )
  
  if (
    file.exists(bigwig_win) &&
    file.info(bigwig_win)$size > 0
  ) {
    
    cat("BigWig already exists — skipping\n")
    
    return(
      tibble(
        GSM = x$GSM,
        sample_id = x$sample_id,
        group = x$group,
        track_type = x$track_type,
        status = "already_exists",
        bigwig_file = bigwig_win,
        bigwig_size_bytes = file.info(bigwig_win)$size,
        error_message = NA_character_
      )
    )
  }
  
  cmd <- sprintf(
    'wsl -d Ubuntu-24.04 -- ~/ucsc_tools/bedGraphToBigWig "%s" ~/ucsc_tools/hg38.chrom.sizes "%s"',
    bedgraph_wsl,
    bigwig_wsl
  )
  
  status <- system(
    cmd
  )
  
  if (
    status == 0 &&
    file.exists(bigwig_win) &&
    file.info(bigwig_win)$size > 0
  ) {
    
    cat(
      "Success |",
      round(file.info(bigwig_win)$size / 1024^2, 3),
      "MB\n"
    )
    
    tibble(
      GSM = x$GSM,
      sample_id = x$sample_id,
      group = x$group,
      track_type = x$track_type,
      status = "success",
      bigwig_file = bigwig_win,
      bigwig_size_bytes = file.info(bigwig_win)$size,
      error_message = NA_character_
    )
    
  } else {
    
    cat("ERROR: conversion failed\n")
    
    tibble(
      GSM = x$GSM,
      sample_id = x$sample_id,
      group = x$group,
      track_type = x$track_type,
      status = "error",
      bigwig_file = bigwig_win,
      bigwig_size_bytes = NA_real_,
      error_message = paste0(
        "bedGraphToBigWig exit status ",
        status
      )
    )
  }
}
test_indices <- c(
  which(regional_manifest$track_type == "WGBS")[1],
  which(regional_manifest$track_type == "ATAC-seq")[1],
  which(regional_manifest$track_type == "H3K27ac")[1]
)

bigwig_test_results <- map_dfr(
  test_indices,
  convert_one_bigwig
)

print(
  as.data.frame(bigwig_test_results)
)

# ============================================================
# Convert all 69 regional bedGraph files to BigWig
# ============================================================

all_bigwig_results <- map_dfr(
  seq_len(nrow(regional_manifest)),
  convert_one_bigwig
)

cat(
  "\nBigWig conversion finished.\n\n"
)

print(
  as.data.frame(
    all_bigwig_results |>
      count(status)
  )
)

total_bigwig_bytes <- sum(
  all_bigwig_results$bigwig_size_bytes,
  na.rm = TRUE
)

cat(
  "\nTotal BigWig size:",
  round(total_bigwig_bytes / 1024^2, 2),
  "MB\n"
)

# ============================================================
# Save BigWig manifest
# ============================================================

bigwig_manifest <- regional_manifest |>
  select(
    GSM,
    sample_id,
    group,
    assay,
    mark,
    track_type,
    title
  ) |>
  left_join(
    all_bigwig_results |>
      select(
        GSM,
        status,
        bigwig_file,
        bigwig_size_bytes
      ),
    by = "GSM"
  ) |>
  arrange(
    track_type,
    group,
    sample_id
  )

write_tsv(
  bigwig_manifest,
  here(
    "metadata",
    "regional_1d_bigwig_tracks.tsv"
  )
)

cat(
  "BigWig manifest saved:\n",
  here("metadata", "regional_1d_bigwig_tracks.tsv"),
  "\n\n"
)

cat(
  "Tracks in BigWig manifest:",
  nrow(bigwig_manifest),
  "\n"
)

cat(
  "Tracks not ready:",
  sum(
    !bigwig_manifest$status %in% c(
      "success",
      "already_exists"
    )
  ),
  "\n"
)