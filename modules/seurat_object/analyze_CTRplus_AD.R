# Script for analyzing combined Seurat object (CTR+ vs AD)

# Load paths from PATHs.txt
source("F:/R_project/PATHs.txt")

# Install and load required packages
required_packages <- c(
  "Seurat", "dplyr", "ggplot2", "patchwork", "sctransform",
  "clustree", "cluster", "glmGamPoi", "presto"
)
install_required_packages(required_packages)

# Create unique paths for current script run
script_name <- "analyze_CTRplus_AD"
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

# Load combined Seurat object
combined_seurat <- readRDS(file.path(r_seurat_merged, "combined_raw_CTRplus_AD.rds"))
log_message(paste("Loaded Seurat object with", dim(combined_seurat)[2], "cells and",
                  dim(combined_seurat)[1], "genes\n"))

# QC metrics visualization pre-filtering
p1 <- VlnPlot(
  combined_seurat, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3, pt.size = 0.1, group.by = "group"
) + ggtitle("QC metrics before filtering")

ggsave(file.path(paths$figures_path, "qc_before_filtering.png"),
       p1, width = 12, height = 6, dpi = 300)

log_message(paste("Saved QC metrics plot before filtering:",
                 file.path(paths$figures_path, "qc_before_filtering.png"), "\n"))

# Filter low quality cells
combined_seurat <- subset(
  combined_seurat,
  subset = nFeature_RNA > 200 &
          nFeature_RNA < 4000 &
          percent.mt < 15
)

log_message(paste("After filtering:", dim(combined_seurat)[2], "cells and",
                 dim(combined_seurat)[1], "genes\n"))

# QC metrics visualization post-filtering
p2 <- VlnPlot(
  combined_seurat, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3, pt.size = 0.1, group.by = "group"
) + ggtitle("QC metrics after filtering")

ggsave(file.path(paths$figures_path, "qc_after_filtering.png"),
       p2, width = 12, height = 6, dpi = 300)

log_message(paste("Saved QC metrics plot after filtering:",
                 file.path(paths$figures_path, "qc_after_filtering.png"), "\n"))

# SCTransform normalization
log_message("Performing SCTransform normalization...\n")
combined_seurat <- SCTransform(
  combined_seurat,
  vars.to.regress = c("percent.mt"),
  verbose = FALSE,
  method = "glmGamPoi"
)
log_message("SCTransform normalization completed\n")

# Save normalized object
norm_obj_path <- file.path(r_seurat_normalized, "combined_normalized_seurat.rds")
temp_obj_path <- file.path(paths$seurat_objects_path, "combined_normalized_seurat.rds")
saveRDS(combined_seurat, norm_obj_path)
saveRDS(combined_seurat, temp_obj_path)
log_message(paste("Normalized object saved:", norm_obj_path, "\n"))
log_message(paste("Temporary copy:", temp_obj_path, "\n"))

# PCA analysis
log_message("Performing PCA analysis...\n")
combined_seurat <- RunPCA(
  combined_seurat,
  features = VariableFeatures(object = combined_seurat),
  npcs = 30
)
log_message("PCA analysis completed\n")

# PCA visualization by groups
p3 <- DimPlot(combined_seurat, reduction = "pca", group.by = "group") +
     ggtitle("PCA separation by groups (CTR+ vs AD)")
ggsave(file.path(paths$figures_path, "pca_groups.png"),
       p3, width = 8, height = 7, dpi = 300)

# PCA visualization by donors
p4 <- DimPlot(combined_seurat, reduction = "pca", group.by = "donor_id") +
     ggtitle("PCA separation by donors")
ggsave(file.path(paths$figures_path, "pca_donors.png"),
       p4, width = 8, height = 7, dpi = 300)

# Elbow plot for optimal PC determination
p5 <- ElbowPlot(combined_seurat, ndims = 30) +
     ggtitle("Elbow plot for optimal PC determination")
ggsave(file.path(paths$figures_path, "pca_elbow.png"),
       p5, width = 8, height = 6, dpi = 300)

# Set optimal PC number for clustering
dims_to_use <- 15  # Typically using 10-15 PCs based on elbow plot and paper
log_message(paste("Using", dims_to_use, "principal components for clustering\n"))

# Build nearest neighbor graph
log_message("Building nearest neighbor graph...\n")
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:dims_to_use)

# Test different resolution values for clustering
log_message("Testing different resolution values for clustering...\n")
resolutions <- seq(0.1, 1.0, by = 0.1)
for (res in resolutions) {
  combined_seurat <- FindClusters(combined_seurat, resolution = res)
}

# Cluster tree visualization
p_clustree <- clustree(combined_seurat, prefix = "SCT_snn_res.")
ggsave(file.path(paths$figures_path, "clustree_resolution.png"),
       p_clustree, width = 10, height = 12, dpi = 300)
log_message(paste("Cluster tree saved:",
                 file.path(paths$figures_path, "clustree_resolution.png"), "\n"))

# Select final resolution
final_resolution <- 0.5
log_message(paste("Selected final resolution:", final_resolution, "\n"))
combined_seurat <- FindClusters(combined_seurat, resolution = final_resolution)

# UMAP visualization
log_message("Running UMAP...\n")
combined_seurat <- RunUMAP(combined_seurat, dims = 1:dims_to_use)
log_message("UMAP completed\n")

# UMAP visualization by clusters
p6 <- DimPlot(combined_seurat, reduction = "umap", label = TRUE, label.size = 5) +
     ggtitle("Microglia clusters")
ggsave(file.path(paths$figures_path, "umap_clusters.png"),
       p6, width = 8, height = 7, dpi = 300)

# UMAP visualization by groups
p7 <- DimPlot(combined_seurat, reduction = "umap", group.by = "group") +
     ggtitle("UMAP separation by groups (CTR+ vs AD)")
ggsave(file.path(paths$figures_path, "umap_groups.png"),
       p7, width = 8, height = 7, dpi = 300)

# UMAP visualization by donors
p8 <- DimPlot(combined_seurat, reduction = "umap", group.by = "donor_id") +
     ggtitle("UMAP separation by donors")
ggsave(file.path(paths$figures_path, "umap_donors.png"),
       p8, width = 8, height = 7, dpi = 300)

# Combined plots of clusters and groups
combined_plot <- p6 + p7
ggsave(file.path(paths$figures_path, "umap_clusters_and_groups.png"),
       combined_plot, width = 16, height = 7, dpi = 300)

# CRITICAL FIX: Prepare object for finding markers after SCTransform
log_message("Preparing SCTransformed data for marker detection...\n")
combined_seurat <- PrepSCTFindMarkers(combined_seurat)

# Find markers for each cluster
log_message("Finding cluster markers...\n")
all_markers <- FindAllMarkers(
  combined_seurat, 
  only.pos = TRUE,
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# Check if markers were found
if (nrow(all_markers) == 0) {
  log_message("WARNING: No cluster markers were found!\n")
} else {
  # Save all markers
  write.csv(all_markers, file.path(paths$results_data_path, "all_cluster_markers.csv"),
            row.names = FALSE)
  log_message(paste("Cluster markers saved:",
                   file.path(paths$results_data_path, "all_cluster_markers.csv"), "\n"))
  
  # Get top 5 markers for each cluster
  log_message(paste("Column names in all_markers:", 
                   paste(colnames(all_markers), collapse=", "), "\n"))
  
  # Use correct column name for cluster information
  if ("cluster" %in% colnames(all_markers)) {
    cluster_col <- "cluster"
  } else {
    cluster_col <- names(all_markers)[1]  # Fallback to first column if cluster not found
  }
  
  log_message(paste("Using column for clustering markers:", cluster_col, "\n"))
  
  top5_markers <- all_markers %>%
                  group_by(!!sym(cluster_col)) %>%
                  top_n(n = 5, wt = avg_log2FC)
  
  # Create heatmap for top markers
  p9 <- DoHeatmap(combined_seurat, features = top5_markers$gene,
                  group.by = "seurat_clusters") +
        ggtitle("Top 5 markers for each cluster")
  ggsave(file.path(paths$figures_path, "top5_markers_heatmap.png"),
         p9, width = 12, height = 10, dpi = 300)
}

# Find differentially expressed genes between CTR+ and AD groups
log_message("Finding differentially expressed genes between CTR+ and AD groups...\n")
# For comparing groups, we can either continue with SCT or switch to RNA assay
# Option 1: Continue with SCT assay (already prepared with PrepSCTFindMarkers)
de_genes <- FindMarkers(
  combined_seurat, 
  ident.1 = "CTR+", 
  ident.2 = "AD",
  group.by = "group", 
  logfc.threshold = 0.1
)

# Check if DE genes were found between groups
if (nrow(de_genes) == 0) {
  log_message("WARNING: No differentially expressed genes found between CTR+ and AD groups with SCT assay.\n")
  log_message("Trying with RNA assay instead...\n")
  
  # Option 2: Switch to RNA assay if SCT comparison fails
  DefaultAssay(combined_seurat) <- "RNA"
  combined_seurat <- NormalizeData(combined_seurat)
  
  de_genes <- FindMarkers(
    combined_seurat, 
    ident.1 = "CTR+", 
    ident.2 = "AD",
    group.by = "group", 
    logfc.threshold = 0.1
  )
  
  # Switch back to SCT assay for the rest of the analysis
  DefaultAssay(combined_seurat) <- "SCT"
}

if (nrow(de_genes) > 0) {
  write.csv(de_genes, file.path(paths$results_data_path, "CTRplus_vs_AD_DE_genes.csv"),
            row.names = TRUE)
  log_message(paste("Differentially expressed genes saved:",
                    file.path(paths$results_data_path, "CTRplus_vs_AD_DE_genes.csv"), "\n"))
  
  # Visualize top 10 differentially expressed genes
  num_de_genes <- min(10, nrow(de_genes))
  top10_de_genes <- rownames(de_genes)[
    order(abs(de_genes$avg_log2FC), decreasing = TRUE)][1:num_de_genes]
  
  p10 <- DotPlot(combined_seurat, features = top10_de_genes, group.by = "group") +
        RotatedAxis() + ggtitle("Top differential genes between CTR+ and AD")
  ggsave(file.path(paths$figures_path, "top10_de_genes_dotplot.png"),
         p10, width = 10, height = 6, dpi = 300)
  
  # Check if specific genes from the paper are in our DE genes list
  paper_genes <- c("SIGLEC1", "CXCL10", "CXCR2")
  paper_de_genes <- paper_genes[paper_genes %in% rownames(de_genes)]
  
  if (length(paper_de_genes) > 0) {
    log_message(paste("Found paper-mentioned genes in DE results:", 
                      paste(paper_de_genes, collapse=", "), "\n"))
    
    p_paper <- VlnPlot(combined_seurat, features = paper_de_genes, 
                       group.by = "group", ncol = length(paper_de_genes)) +
      ggtitle("DE genes mentioned in paper")
    ggsave(file.path(paths$figures_path, "paper_de_genes.png"),
           p_paper, width = 4 * length(paper_de_genes), height = 6, dpi = 300)
  }
}

# Known microglia and neurodegeneration markers
known_markers <- c(
  "P2RY12", "TMEM119", "CX3CR1", "ITGAX", "AXL", "CLEC7A",
  "APOE", "TREM2", "CD68", "CST7", "CD74", "HLA-DRA"
)

found_markers <- known_markers[known_markers %in% rownames(combined_seurat)]
if (length(found_markers) > 0) {
  p11 <- FeaturePlot(
    combined_seurat,
    features = found_markers,
    ncol = 4, pt.size = 0.1
  )
  ggsave(file.path(paths$figures_path, "known_markers_featureplot.png"),
         p11, width = 16, height = 12, dpi = 300)
}

# Save final Seurat object
final_obj_path <- file.path(r_seurat_merged, "combined_final_CTRplus_AD.rds")
temp_obj_path <- file.path(paths$seurat_objects_path, "combined_final_CTRplus_AD.rds")
saveRDS(combined_seurat, final_obj_path)
saveRDS(combined_seurat, temp_obj_path)
log_message(paste("Final object saved:", final_obj_path, "\n"))
log_message(paste("Temporary copy:", temp_obj_path, "\n"))

log_message(paste("\nAnalysis completed:", Sys.time(), "\n"))
log_message(paste("All results saved to:", paths$results_path, "\n"))