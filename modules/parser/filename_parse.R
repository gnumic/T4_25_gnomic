# Script for analyzing file name structure in directory

# Load paths from PATHs.txt
source("F:/R_project/PATHs.txt")

# Install and load required packages
required_packages <- c("data.table")
install_required_packages(required_packages)

# Create unique paths for current script run
script_name <- "filename_parse"
paths <- create_unique_output_paths(script_name)

# Start logging
log_file <- paths$log_file
cat(paste("Script execution start:", Sys.time(), "\n"), file = log_file)
cat(paste("Results will be saved to:", paths$results_path, "\n"), 
    file = log_file, append = TRUE)

# Function to duplicate output to console and log
log_message <- function(message) {
  cat(message)
  cat(message, file = log_file, append = TRUE)
}

# Get list of all files in directory
all_files <- list.files(r_raw_data, pattern = "\\.txt$", full.names = TRUE)
log_message(paste("Found", length(all_files), "text files in data directory\n"))

# Create empty data.frame for results
file_components <- data.frame(
  filename = character(),
  full_path = character(),
  GSM_ID = character(),
  donor_ID = character(),
  condition = character(),
  sample_ID = character(),
  seq_ID = character(),
  stringsAsFactors = FALSE
)

# Process all files
for (file_path in all_files) {
  filename <- basename(file_path)
  base_name <- sub("_completeCounts\\.txt$", "", filename)
  parts <- strsplit(base_name, "_")[[1]]
  
  GSM_ID <- parts[1]
  donor_ID <- ifelse(length(parts) > 1, parts[2], NA)
  condition <- ifelse(length(parts) > 2, parts[3], NA)
  sample_ID <- ifelse(length(parts) > 3, parts[4], NA)
  seq_ID <- ifelse(length(parts) > 4, parts[5], NA)
  
  file_components <- rbind(file_components, data.frame(
    filename = filename,
    full_path = file_path,
    GSM_ID = GSM_ID,
    donor_ID = donor_ID,
    condition = condition,
    sample_ID = sample_ID,
    seq_ID = seq_ID,
    stringsAsFactors = FALSE
  ))
}

# Count files by condition
conditions_summary <- table(file_components$condition)
log_message("\nNumber of files by condition:\n")
print(conditions_summary)
capture.output(conditions_summary, file = log_file, append = TRUE)

# Count files by donor
donors_summary <- table(file_components$donor_ID)
log_message("\nNumber of files by donor:\n")
print(donors_summary)
capture.output(donors_summary, file = log_file, append = TRUE)

# Save results to CSV file
csv_path <- file.path(paths$results_data_path, "file_components.csv")
write.csv(file_components, csv_path, row.names = FALSE)
log_message(paste("\nFile components saved to:", csv_path, "\n"))

# Save results to RDS file
rds_path <- file.path(paths$results_data_path, "file_components.rds")
saveRDS(file_components, rds_path)
log_message(paste("File components saved to:", rds_path, "\n"))

# Create list with name components
file_parts_list <- list()
for (i in seq_len(nrow(file_components))) {
  file_name <- file_components$filename[i]
  file_parts_list[[file_name]] <- list(
    GSM_ID = file_components$GSM_ID[i],
    donor_ID = file_components$donor_ID[i],
    condition = file_components$condition[i],
    sample_ID = file_components$sample_ID[i],
    seq_ID = file_components$seq_ID[i]
  )
}

# Save list for future use
parts_list_path <- file.path(paths$results_data_path, "file_parts_list.rds")
saveRDS(file_parts_list, parts_list_path)
log_message(paste("File parts list saved to:", parts_list_path, "\n"))

# Extract files for specific donor (example)
donor_files <- file_components[file_components$donor_ID == "2018-018", ]
log_message(paste("\nFound", nrow(donor_files), "files for donor 2018-018\n"))

# Extract CTR and AD samples
ctr_files <- file_components[file_components$condition == "CTR", ]
log_message(paste("Found", nrow(ctr_files), "Control (CTR) samples\n"))

ad_files <- file_components[file_components$condition == "AD", ]
log_message(paste("Found", nrow(ad_files), "Alzheimer's (AD) samples\n"))

log_message(paste("\nScript execution completed:", Sys.time(), "\n"))