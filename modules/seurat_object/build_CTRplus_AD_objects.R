# Script for creating Seurat objects for CTR+ and AD groups using parsed file information

# Load paths from PATHs.txt
source("F:/R_project/PATHs.txt")

# Install and load required packages
required_packages <- c(
  "Seurat", "dplyr", "ggplot2", "data.table",
  "sctransform", "glmGamPoi"
)
install_required_packages(required_packages)

# Create unique paths for current script run
script_name <- "build_CTRplus_AD_objects"
paths <- create_unique_output_paths(script_name)

# Create folders for storing Seurat objects if they don't exist
dir.create(r_seurat_raw, recursive = TRUE, showWarnings = FALSE)
dir.create(r_seurat_merged, recursive = TRUE, showWarnings = FALSE)

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

# Function to find the most recent file_components.rds
find_latest_parsed_data <- function() {
  # List all directories matching the pattern
  all_parse_dirs <- list.dirs(path = r_export, recursive = FALSE)
  parse_dirs <- all_parse_dirs[grepl("filename_parse_\\d+_\\d+$", all_parse_dirs)]
  
  if (length(parse_dirs) == 0) {
    stop("No filename_parse results found. Run filename_parse.R first.")
  }
  
  # Get the most recent directory
  latest_dir <- sort(parse_dirs, decreasing = TRUE)[1]
  
  # Return path to the file_components.rds
  file.path(latest_dir, "results", "file_components.rds")
}

# 1. Load parsed file information
log_message("Loading parsed file information...\n")
parsed_data_file <- find_latest_parsed_data()
file_components <- readRDS(parsed_data_file)
log_message(paste("Loaded file components from:", parsed_data_file, "\n"))
log_message(paste("Total files parsed:", nrow(file_components), "\n"))

# 2. Read donors file to determine CTR+ and AD groups
donors_info <- read.csv(
  file.path(r_supplementary, "donors.csv"),
  sep = ";", skip = 1, header = TRUE
)

# Select donors from CTR+ and AD groups
ctrplus_donors <- donors_info$DonorID[donors_info$Group == "CTR +"]
ad_donors <- donors_info$DonorID[donors_info$Group == "AD"]

log_message(paste("CTR+ group donors:", paste(ctrplus_donors, collapse = ", "), "\n"))
log_message(paste("AD group donors:", paste(ad_donors, collapse = ", "), "\n"))

# 3. Read metadata file
metadata <- read.delim(
  file.path(r_my_data, "SingleCellsAlzheimers 2025-04-16 09.15.tsv"),
  sep = "\t", header = TRUE, quote = ""
)

# 4. Function to create unique gene names
make_unique_names <- function(names) {
  counts <- table(names)
  duplicates <- names(counts[counts > 1])
  
  for (dup in duplicates) {
    indices <- which(names == dup)
    for (i in seq_along(indices)) {
      if (i > 1) {
        names[indices[i]] <- paste0(dup, ".", i)
      }
    }
  }
  names
}

# 5. Function to create Seurat object from expression file
create_seurat_from_file <- function(file_path, donor_id, group) {
  file_name <- basename(file_path)
  
  log_message(paste("Processing file:", file_name, "for donor:", donor_id, "group:", group, "\n"))
  
  # Read expression data
  counts <- data.table::fread(file_path, header = TRUE)
  
  # Convert to matrix
  gene_ids <- counts[[1]]
  gene_names <- make_unique_names(gene_ids)
  counts_matrix <- as.matrix(counts[, -1])
  rownames(counts_matrix) <- gene_names
  
  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts_matrix,
    project = donor_id,
    min.cells = 3,
    min.features = 200
  )
  
  # Add metadata
  seurat_obj$donor_id <- donor_id
  seurat_obj$group <- group
  seurat_obj$file_name <- file_name
  
  # Calculate percentage of mitochondrial genes
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    seurat_obj, pattern = "^MT-"
  )
  
  # Add additional metadata from file
  donor_metadata <- metadata[grep(
    donor_id, metadata$donor_organism.biomaterial_core.biomaterial_id
  ), ]
  
  if (nrow(donor_metadata) > 0) {
    seurat_obj$sex <- donor_metadata$donor_organism.sex[1]
    seurat_obj$age <- donor_metadata$donor_organism.organism_age[1]
    seurat_obj$disease <- donor_metadata$donor_organism.diseases[1]
  }
  
  seurat_obj
}

# 6. Create Seurat objects for each group using parsed file information
ctrplus_objects <- list()
ad_objects <- list()

# Filter files for CTR+ donors
log_message("\nProcessing files for CTR+ group\n")
for (donor in ctrplus_donors) {
  # Use parsed data to filter files
  donor_files_df <- file_components[file_components$donor_ID == donor & 
                                   file_components$condition == "CTR", ]
  
  if (nrow(donor_files_df) > 0) {
    for (i in seq_len(nrow(donor_files_df))) {
      file_path <- donor_files_df$full_path[i]
      seurat_obj <- create_seurat_from_file(file_path, donor, "CTR+")
      
      if (!is.null(seurat_obj)) {
        # Extract sample_ID from the parsed data
        sample_id <- donor_files_df$sample_ID[i]
        obj_name <- paste0(donor, "_", sample_id)
        ctrplus_objects[[obj_name]] <- seurat_obj
        
        # Save individual object in both directories
        raw_obj_path <- file.path(r_seurat_raw, paste0("CTRplus_", obj_name, ".rds"))
        temp_obj_path <- file.path(paths$seurat_objects_path, paste0("CTRplus_", obj_name, ".rds"))
        saveRDS(seurat_obj, raw_obj_path)
        saveRDS(seurat_obj, temp_obj_path)
        
        log_message(paste("Object saved:", raw_obj_path, "\n"))
      }
    }
  }
}

# Filter files for AD donors
log_message("\nProcessing files for AD group\n")
for (donor in ad_donors) {
  # Use parsed data to filter files
  donor_files_df <- file_components[file_components$donor_ID == donor & 
                                   file_components$condition == "AD", ]
  
  if (nrow(donor_files_df) > 0) {
    for (i in seq_len(nrow(donor_files_df))) {
      file_path <- donor_files_df$full_path[i]
      seurat_obj <- create_seurat_from_file(file_path, donor, "AD")
      
      if (!is.null(seurat_obj)) {
        # Extract sample_ID from the parsed data
        sample_id <- donor_files_df$sample_ID[i]
        obj_name <- paste0(donor, "_", sample_id)
        ad_objects[[obj_name]] <- seurat_obj
        
        # Save individual object in both directories
        raw_obj_path <- file.path(r_seurat_raw, paste0("AD_", obj_name, ".rds"))
        temp_obj_path <- file.path(paths$seurat_objects_path, paste0("AD_", obj_name, ".rds"))
        saveRDS(seurat_obj, raw_obj_path)
        saveRDS(seurat_obj, temp_obj_path)
        
        log_message(paste("Object saved:", raw_obj_path, "\n"))
      }
    }
  }
}

log_message(paste("\nCreated", length(ctrplus_objects), "Seurat objects for CTR+ group\n"))
log_message(paste("Created", length(ad_objects), "Seurat objects for AD group\n"))

# 7. Merge CTR+ objects into one Seurat object
if (length(ctrplus_objects) > 0) {
  if (length(ctrplus_objects) == 1) {
    combined_ctrplus <- ctrplus_objects[[1]]
  } else {
    combined_ctrplus <- merge(
      ctrplus_objects[[1]],
      y = ctrplus_objects[2:length(ctrplus_objects)],
      add.cell.ids = names(ctrplus_objects),
      project = "CTRplus"
    )
  }
  
  # Save merged CTR+ object
  merged_obj_path <- file.path(r_seurat_merged, "combined_CTRplus.rds")
  temp_obj_path <- file.path(paths$seurat_objects_path, "combined_CTRplus.rds")
  saveRDS(combined_ctrplus, merged_obj_path)
  saveRDS(combined_ctrplus, temp_obj_path)
  
  log_message(paste("\nSaved combined CTR+ object with", ncol(combined_ctrplus), "cells\n"))
  log_message(paste("Path:", merged_obj_path, "\n"))
}

# 8. Merge AD objects into one Seurat object
if (length(ad_objects) > 0) {
  if (length(ad_objects) == 1) {
    combined_ad <- ad_objects[[1]]
  } else {
    combined_ad <- merge(
      ad_objects[[1]],
      y = ad_objects[2:length(ad_objects)],
      add.cell.ids = names(ad_objects),
      project = "AD"
    )
  }
  
  # Save merged AD object
  merged_obj_path <- file.path(r_seurat_merged, "combined_AD.rds")
  temp_obj_path <- file.path(paths$seurat_objects_path, "combined_AD.rds")
  saveRDS(combined_ad, merged_obj_path)
  saveRDS(combined_ad, temp_obj_path)
  
  log_message(paste("\nSaved combined AD object with", ncol(combined_ad), "cells\n"))
  log_message(paste("Path:", merged_obj_path, "\n"))
}

# 9. Merge CTR+ and AD objects
if (exists("combined_ctrplus") && exists("combined_ad")) {
  combined_all <- merge(
    combined_ctrplus, 
    y = combined_ad,
    add.cell.ids = c("CTRplus", "AD"),
    project = "CTRplus_vs_AD"
  )
  
  # Save combined object
  merged_obj_path <- file.path(r_seurat_merged, "combined_raw_CTRplus_AD.rds")
  temp_obj_path <- file.path(paths$seurat_objects_path, "combined_raw_CTRplus_AD.rds")
  saveRDS(combined_all, merged_obj_path)
  saveRDS(combined_all, temp_obj_path)
  
  log_message(paste("\nSaved combined object with", ncol(combined_all), "cells\n"))
  log_message(paste("Path:", merged_obj_path, "\n"))
}

log_message(paste("\nSeurat object creation process completed:", Sys.time(), "\n"))
log_message(paste("All results saved to:", paths$results_path, "\n"))