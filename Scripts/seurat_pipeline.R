# Author: Joseph Aston Phiri
# Description: Modular Seurat single-cell RNA-seq pipeline for multiple samples
# Version: 1.0
# Date: 2025-01-10

# ------------------------------- #
#        Required Libraries       #
# ------------------------------- #
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(harmony)

# ------------------------------- #
#       Process Single Sample     #
# ------------------------------- #
process_single_sample <- function(sample_name, data_dir = "Data/Single_Cell_Data",
                                  min_features = 30,
                                  filter_params = list(nUMI = c(100, 50000), nGene = 100, mitoRatio = 0.30)) {
  
  # Load and create Seurat object
  raw_data <- Read10X(data.dir = file.path(data_dir, paste0(sample_name, "_raw_feature_bc_matrix")))
  seurat_obj <- CreateSeuratObject(counts = raw_data, min.features = min_features, project = sample_name)
  
  # Compute QC metrics
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  seurat_obj$mitoRatio <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") / 100
  
  # Rename metadata
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    dplyr::rename(seq_folder = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  
  # Filtering
  seurat_filt <- subset(seurat_obj,
                        subset = nUMI >= filter_params$nUMI[1] &
                          nUMI <= filter_params$nUMI[2] &
                          nGene >= filter_params$nGene &
                          mitoRatio < filter_params$mitoRatio)
  
  # Standard workflow
  seurat_filt <- NormalizeData(seurat_filt)
  seurat_filt <- FindVariableFeatures(seurat_filt)
  seurat_filt <- ScaleData(seurat_filt)
  seurat_filt <- RunPCA(seurat_filt)
  seurat_filt <- FindNeighbors(seurat_filt)
  seurat_filt <- FindClusters(seurat_filt)
  seurat_filt <- RunUMAP(seurat_filt, dims = 1:30)
  
  return(seurat_filt)
}

# ------------------------------- #
#     Run Pipeline on All Files   #
# ------------------------------- #
run_all_samples <- function(sample_list, mito_dict) {
  sample_objs <- list()
  for (sample in sample_list) {
    cat("Processing", sample, "\n")
    mito_cutoff <- mito_dict[[sample]]
    sample_objs[[sample]] <- process_single_sample(sample_name = sample,
                                                   filter_params = list(nUMI = c(100, 50000),
                                                                        nGene = 100,
                                                                        mitoRatio = mito_cutoff))
  }
  return(sample_objs)
}

# ------------------------------- #
#     Merge All Sample Objects    #
# ------------------------------- #
merge_all_samples <- function(sample_objs, save_path = "Data/Single_Cell_Data/all_merged.RData") {
  merged <- Reduce(function(x, y) merge(x, y, add.cell.ids = NULL), sample_objs)
  save(merged, file = save_path)
  return(merged)
}

# ------------------------------- #
#       Harmony Integration       #
# ------------------------------- #
run_harmony_integration <- function(seurat_obj, group.by = "seq_folder", dims = 1:30) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 50)
  
  seurat_obj <- IntegrateLayers(object = seurat_obj,
                                method = HarmonyIntegration,
                                group.by = group.by,
                                dims = dims,
                                orig.reduction = "pca",
                                new.reduction = "harmony",
                                verbose = FALSE)
  
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.3, cluster.name = "harmony_clusters", method = "igraph")
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = dims, reduction.name = "umap.harmony")
  seurat_obj <- RunTSNE(seurat_obj, reduction = "harmony", dims = dims, reduction.name = "tsne.harmony")
  
  return(seurat_obj)
}

# ------------------------------- #
#     Example Usage (Pipeline)    #
# ------------------------------- #
# Define sample names and mitoRatio cutoffs
sample_list <- c("CUH124", "CUG11X", "CUF131", "CUF130", "CUF134", "CUF13I", "CUF135",
                 "CUF136", "CUF13J", "CUF13K", "CUF137")

mito_dict <- list(
  CUH124 = 0.25, CUG11X = 0.20, CUF130 = 0.30,
  CUF131 = 0.30, CUF134 = 0.25, CUF135 = 0.25,
  CUF136 = 0.30, CUF13J = 0.30, CUF13K = 0.30,
  CUF137 = 0.20, CUF13I = 0.30
)

# Run and merge
all_samples <- run_all_samples(sample_list, mito_dict)
all_merged <- merge_all_samples(all_samples)

# Run Harmony integration
all_merged <- run_harmony_integration(all_merged)

# Save integrated object
save(all_merged, file = "Data/Single_Cell_Data/all_merged_harmony.RData")

# Plot UMAP
DimPlot(all_merged, reduction = "umap.harmony", label = TRUE, repel = TRUE) + theme_bw()

