library(data.table)
library(stringr)
library(optparse)
set.seed(42)

# Define command line arguments
option_list <- list(
  make_option("--data", type="character", 
              help="Path to file containing list of data files"),
  make_option("--output_dir", type="character", default=".",
              help="Output directory for split files"),
  make_option("--min_num_seqs", type="integer", default=20,
              help="Minimum number of sequences required for an organism [default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$data)) {
  stop("Usage: Rscript script.R --data <path_to_data_files_list> [--output_dir <output_directory>]")
}

input_file <- opt$data
output_dir <- opt$output_dir
min_num_seqs <- opt$min_num_seqs
# Validate input file exists
if (!file.exists(input_file)) {
  stop(sprintf("Input file '%s' does not exist", input_file))
}

# Read and clean file paths
data_files <- str_trim(readLines(input_file))
data_files <- data_files[str_length(data_files) > 0]  # Remove empty lines

# Validate we have files to process
if (length(data_files) == 0) {
  stop("No valid file paths found in input file")
}

message("Found ", length(data_files), " files to process")

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

find_chromosomes_for_target <- function(dt_chr_counts, target_size, excluded_chrs = character()) {
  dt <- copy(dt_chr_counts[!chr %in% excluded_chrs])
  selected_chrs <- character()
  current_count <- 0
  
  # Sort by count to try smaller chromosomes first
  setorder(dt, N)
  
  while (current_count < target_size && nrow(dt) > 0) {
    # Try to find a combination that gets us closer to exactly 5%
    remaining_target <- target_size - current_count
    
    # Find the chromosome that gets us closest to target without going under 5%
    dt[, diff := N - remaining_target]
    best_match <- dt[diff >= 0][which.min(diff)]
    
    # If no chromosome gets us above target with a single addition,
    # take the smallest remaining chromosome
    if (nrow(best_match) == 0) {
      best_match <- dt[1]
    }
    
    selected_chrs <- c(selected_chrs, best_match$chr)
    current_count <- current_count + best_match$N
    dt <- dt[chr != best_match$chr]
    
    # If we're already above 5%, check if removing the last added chromosome
    # would get us closer to 5%
    if (length(selected_chrs) > 1 && current_count > target_size) {
      test_count <- current_count - best_match$N
      if (test_count >= target_size && 
          abs(test_count - target_size) < abs(current_count - target_size)) {
        selected_chrs <- selected_chrs[-length(selected_chrs)]
        current_count <- test_count
      }
    }
  }
  
  list(selected_chrs = selected_chrs, count = current_count)
}

get_chr_splits <- function(dt) {
  chr_counts <- dt[, .(N = .N), by = chr]
  total_seqs <- sum(chr_counts$N)
  target_size <- total_seqs * 0.05  # Exactly 5%
  
  # Get validation chromosomes
  val_result <- find_chromosomes_for_target(chr_counts, target_size)
  
  # Get test chromosomes
  test_result <- find_chromosomes_for_target(chr_counts, target_size, val_result$selected_chrs)
  
  # Calculate percentages
  list(
    val_chr = val_result$selected_chrs,
    test_chr = test_result$selected_chrs,
    val_pct = (val_result$count / total_seqs) * 100,
    test_pct = (test_result$count / total_seqs) * 100,
    train_pct = ((total_seqs - val_result$count - test_result$count) / total_seqs) * 100
  )
}

get_stats <- function(dt) {
  num_genes = dt[, gene_id %>% unique() %>% length()]
  num_trs = dt[, transcript_id %>% unique() %>% length()]
  num_tokens = dt[, sum(seq_len)] %>% sum()
  num_exonic_tokens = dt[, str_count(text, '[ACGU]')] %>% sum()
  num_intronic_tokens = dt[, str_count(text, '[acgu]')] %>% sum()
  num_unk_tokens = dt[, str_count(text, '[Nn]')] %>% sum()
  num_contains_unk = dt[, sum(contains_unk)]
  
  list(
    num_genes = num_genes,
    num_trs = num_trs,
    num_tokens = num_tokens,
    num_exonic_tokens = num_exonic_tokens,
    num_intronic_tokens = num_intronic_tokens,
    num_unk_tokens = num_unk_tokens,
    num_contains_unk = num_contains_unk
  )
}

# Initialize metadata list
metadata_list <- vector("list", length(data_files))

# Process each file
for (i in seq_along(data_files)) {
  data_file <- data_files[i]
  message(str_c(rep("-", 30), collapse = ""))
  message("Processing ", data_file, "...")
  
  # Extract organism name
  organism <- str_replace(basename(data_file), "\\.csv\\.gz$", "")

  # Save splits
  output_files <- str_c(
    file.path(output_dir, organism),
    c(".train.csv.gz", ".validation.csv.gz", ".test.csv.gz")
  )
  
  # Check if output files already exist
  if (all(file.exists(output_files))) {
    message("Output files already exist for ", organism, ", skipping...")
    next
  }

  # Read data
  dt <- fread(data_file)
  # Skip organisms with too few sequences
  if (nrow(dt) < min_num_seqs) {
    message(sprintf("Skipping %s - only %d sequences (minimum %d required)", 
                   organism, nrow(dt), min_num_seqs))
    next
  }

  if ("chr" %in% names(dt)) {
    # Get splits
    splits <- get_chr_splits(dt)
    
    # Create split datasets
    dt_val <- dt[chr %in% splits$val_chr]
    dt_test <- dt[chr %in% splits$test_chr]
    dt_train <- dt[!chr %in% c(splits$val_chr, splits$test_chr)]
    
    # Calculate stats for each split
    train_stats <- get_stats(dt_train)
    val_stats <- get_stats(dt_val)
    test_stats <- get_stats(dt_test)
    
    # Write files with compression
    fwrite(dt_train, output_files[1], compress = "gzip")
    fwrite(dt_val, output_files[2], compress = "gzip")
    fwrite(dt_test, output_files[3], compress = "gzip")
    
    # Store metadata with stats
    metadata_splits <- list(
      list(
        output_file = output_files[1],
        organism = organism,
        split = "train",
        chr = "remaining",
        pct = splits$train_pct
      ),
      list(
        output_file = output_files[2],
        organism = organism,
        split = "validation",
        chr = str_c(splits$val_chr, collapse = ","),
        pct = splits$val_pct
      ),
      list(
        output_file = output_files[3],
        organism = organism,
        split = "test",
        chr = str_c(splits$test_chr, collapse = ","),
        pct = splits$test_pct
      )
    )
    
    # Add stats to each split
    metadata_splits[[1]] <- c(metadata_splits[[1]], train_stats)
    metadata_splits[[2]] <- c(metadata_splits[[2]], val_stats)
    metadata_splits[[3]] <- c(metadata_splits[[3]], test_stats)
    
    # Convert to data.table and store
    metadata_list[[i]] <- rbindlist(lapply(metadata_splits, as.data.table))
    
    # Print progress
    message(sprintf(
      "Split sizes - Train: %.2f%%, Validation: %.2f%%, Test: %.2f%%",
      splits$train_pct, splits$val_pct, splits$test_pct
    ))
  }
}

# Combine and save metadata
metadata <- rbindlist(metadata_list)
fwrite(metadata, file.path(output_dir, "split_metadata.csv"), append = TRUE)
