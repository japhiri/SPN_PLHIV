
# scRNA-seq Modular Seurat Pipeline
# Author: Joseph Aston Phiri
# Description: Modular pipeline for processing, QC, filtering, and merging scRNA-seq samples using Seurat

library(Seurat)
library(dplyr)
library(stringr)
library(patchwork)
library(harmony)
library(presto)
library(BPCells)
library(tidyverse)
library(SingleR)
library(celldex)

# Function to process a single sample
process_single_sample <- function(sample_name, data_dir = "Data/", 
                                  min_features = 30,
                                  filter_params = list(nCount_RNA = c(100, 50000), nFeature_RNA = 100, mitoRatio = 0.30)) {

  seurat_data <- Read10X(data.dir = file.path(data_dir, sample_name))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.features = min_features, project = sample_name)

  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-") / 100

  filt_obj <- subset(seurat_obj,
                     subset = nCount_RNA >= filter_params$nCount_RNA[1] & 
                       nCount_RNA <= filter_params$nCount_RNA[2] & 
                       nFeature_RNA >= filter_params$nFeature_RNA & 
                              mitoRatio < filter_params$mitoRatio)

  return(filt_obj)
}

# Define sample list and filtering parameters
sample_list <- c("CUH124", "CUG11X", "CUF131", "CUF130", "CUF134", "CUF13I", "CUF135", 
                 "CUF136", "CUF13J", "CUF13K", "CUF137", "CUF12Y"
                 )

mito_dict <- list(
  CUH124 = 0.25, CUG11X = 0.20, CUF130 = 0.30,
  CUF131 = 0.30, CUF134 = 0.25, CUF135 = 0.25,
  CUF136 = 0.30, CUF13J = 0.30, CUF13K = 0.30,
  CUF137 = 0.20, CUF12Y = 0.30, CUF13I = 0.30
)

sample_objs <- list()

for (sample in sample_list) {
  mito_cutoff <- mito_dict[[sample]]
  sample_name_raw <- paste0(sample, "_raw_feature_bc_matrix")
  sample_objs[[sample]] <- process_single_sample(sample_name = sample_name_raw,
                                                 filter_params = list(nCount_RNA = c(100, 50000),
                                                                      nFeature_RNA = 100,
                                                                      mitoRatio = mito_cutoff))
}

# Merge all samples
all_merged <- Reduce(function(x, y) merge(x, y, add.cell.ids = NULL), sample_objs)
# Remove Sample CUF12Y due to poor quality
all_merged <- subset(all_merged,cells = which(all_merged$orig.ident!="CUF12Y_raw_feature_bc_matrix"))

# Create metadata columns
annotate_seurat_metadata <- function(seurat_obj) {
  metadata <- seurat_obj@meta.data
  
  # Define mapping
  mapping <- list(
    sample = c(
      "CUH124", "CUG11X", "CUF131", "CUF130", "CUF134", 
      "CUF13I", "CUF135", "CUF136", "CUF13J", "CUF13K", "CUF137"
    ),
    HIV_Status = c(
      "CUH124" = "PLHIV-ART>1yr", "CUG11X" = "PLHIV-ART>1yr", "CUF131" = "HIV-",
      "CUF130" = "PLHIV-ART>1yr", "CUF134" = "PLHIV-ART<3m", "CUF13I" = "PLHIV-ART>1yr",
      "CUF135" = "PLHIV-ART>1yr", "CUF136" = "HIV-", "CUF13J" = "PLHIV-ART<3m",
      "CUF13K" = "PLHIV-ART<3m", "CUF137" = "HIV-"
    ),
    Carriage_Status = c(
      "CUH124" = "SPN+", "CUG11X" = "SPN-", "CUF131" = "SPN+", "CUF130" = "SPN-",
      "CUF134" = "SPN+", "CUF13I" = "SPN+", "CUF135" = "SPN+", "CUF136" = "SPN+",
      "CUF13J" = "SPN+", "CUF13K" = "SPN-", "CUF137" = "SPN+"
    ),
    Carriage_density = c(
      "CUH124" = 5862500, "CUF131" = 26800000, "CUF134" = 686750,
      "CUF13I" = 670, "CUF135" = 33500, "CUF136" = 67, "CUF13J" = 26800,
      "CUF137" = 10720
    ),
    serotype = c(
      "CUH124" = "NVT", "CUF131" = "3", "CUF134" = "9", "CUF13I" = "1",
      "CUF135" = "10", "CUF136" = "NVT", "CUF13J" = "10", "CUF137" = "NVT"
    ),
    age = c(
      "CUH124" = 35, "CUG11X" = 21, "CUF131" = 28, "CUF130" = 34, "CUF134" = 26,
      "CUF13I" = 27, "CUF135" = 34, "CUF136" = 36, "CUF13J" = 24,
      "CUF13K" = 27, "CUF137" = 28
    ),
    sex = c(
      "CUH124" = "Female", "CUG11X" = "Female", "CUF131" = "Female", "CUF130" = "Female",
      "CUF134" = "Male", "CUF13I" = "Male", "CUF135" = "Male", "CUF136" = "Female",
      "CUF13J" = "Female", "CUF13K" = "Female", "CUF137" = "Male"
    ),
    Viral_Load = c(
      "CUH124" = "Not Detected", "CUG11X" = "Not Detected", "CUF131" = "Not Detected",
      "CUF130" = "Not Detected", "CUF134" = 312500, "CUF13I" = "Not Detected",
      "CUF135" = "Not Detected", "CUF136" = "Not Detected", "CUF13J" = "Not Detected",
      "CUF13K" = "Not Detected", "CUF137" = "Not Detected"
    ),
    abs_CD4_count = c(
      "CUH124" = 805, "CUG11X" = 488, "CUF131" = 381, "CUF130" = 535,
      "CUF134" = 14, "CUF13I" = 452, "CUF135" = 274, "CUF136" = 1179,
      "CUF13J" = 689, "CUF13K" = 907, "CUF137" = 431
    )
  )
  
  # Add new columns to metadata
  for (col in names(mapping)) {
    metadata[[col]] <- NA
    map <- mapping[[col]]
    if (is.null(names(map))) {
      for (id in map) {
        metadata[[col]][str_detect(metadata$orig.ident, paste0("^", id))] <- id
      }
    } else {
      for (id in names(map)) {
        metadata[[col]][str_detect(metadata$orig.ident, paste0("^", id))] <- map[[id]]
      }
    }
  }
  
  # Assign updated metadata back to Seurat object
  seurat_obj@meta.data <- metadata
  return(seurat_obj)
}
# Annotate Seurat object
all_merged <- annotate_seurat_metadata(all_merged)

# View updated metadata
head(all_merged@meta.data)

# Standard seurat pipeline
all_merged <- NormalizeData(all_merged)
all_merged <- FindVariableFeatures(all_merged)
all_merged <- ScaleData(all_merged)
#all_merged <- SCTransform(all_merged,  vars.to.regress = c("nCount_RNA","mitoRatio","log10GenesPerUMI"), variable.features.n = 3000, verbose = TRUE,return.only.var.genes = FALSE)
all_merged <- RunPCA(all_merged)

# Batch correction
all_merged <- harmony::RunHarmony(all_merged, group.by.vars = c("orig.ident"), reduction = "pca", 
                                  dims.use = 1:30, assay.use = "RNA", reduction.save = "harmony")
all_merged<-FindNeighbors(all_merged, reduction = "harmony", dims = 1:30, assay = "RNA")
all_merged<-FindClusters(all_merged, resolution = 0.3, cluster.name = "harmony.clusters")
set.seed(1234)
all_merged<-RunUMAP(all_merged, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", assay = "RNA", set.seed = 1234))
# Dimplot
DimPlot(all_merged, reduction = "umap.harmony", label = TRUE)
FeaturePlot(all_merged, features = c("MUC5AC"), reduction = "umap.harmony")
save(all_merged, file = "HIV-PAPER/outputs/all_merged.RData")
load("HIV-PAPER/outputs/all_merged.RData")
# Find Markers for each cluster
all_merged <- JoinLayers(all_merged)
all_merged_markers <- FindAllMarkers(all_merged, 
                                     assay = "RNA",
                                     logfc.threshold = 0.25,
                                     test.use = "wilcox",
                                     slot = "data",
                                     min.pct = 0.1,
                                     min.diff.pct = -Inf,
                                     verbose = TRUE,
                                     only.pos = TRUE,
                                     max.cells.per.ident = Inf,
                                     random.seed = 1,
                                     min.cells.feature = 3,
                                     min.cells.group = 3,
                                     mean.fxn = NULL,
                                     fc.name = NULL,
                                     base = 2,
                                     return.thresh = 0.01,
                                     densify = FALSE)
write_csv(all_merged_markers,"HIV-PAPER/Results/all_merged_markers.csv")
# RenameIdents
all_merged <- RenameIdents(all_merged,
                           "0" = "0",
                           "1" = "1",
                           "2" = "Basal cells",
                           "3" = "Secretory cells",
                           "4" = "Goblet cells",
                           "5" = "T cells",
                           "6" = "FOXJ1++ Ciliated cells",
                           "7" = "Phagocytes",
                           "8" = "Ciliated cells",
                           "9" = "Ciliated cells",
                           "10" = "Neutrophils",
                           "11" = "Squamous cells",
                           "12" = "Deuterosomal cells",
                           "13" = "B cells",
                           "14" = "Ciliated cells",
                           "15" = "Ionocytes")
# Dimplot
DimPlot(all_merged, 
        reduction = "umap.harmony", 
        label = TRUE)+
  NoLegend()

VlnPlot(all_merged, 
        features = c("KRT5","TP63","ITGA6",
                     "MUC2","TFF3","AGR2",
                     "SPDEF","SCGB1A1","BPIFB1",
                     "SPINK1"),
        ncol = 3)

celltype_markers <- c(
  "CYP2F1","SERPINB3","MUC5AC",
  "SPDEF","VMO1","AQP5","BPIFB1",
  "PTPRC","CD3E","CD3D","CD8A","CD4",
  "TYROBP","HLA-DPA1","HLA-DPB1",
  "MS4A1","CD79A","CD19",
  "CAPS","CFAP157","DNAAF1","FOXJ1",
  "G0S2","CXCL8","CSF3R","KRT4",
  "SPRR3","SPRR2A","KRT6A",
  "CDC20B","CCNO","MSMB",
  "CFTR","SCNN1B","RARRES2",
  "CD83","CD86","LYZ","FCER1G",
  "PCP4L1","CSRP2","SERPINF1","ITGA6",
  "KRT5","CD1C","ITGAX","XCR1","CD68","CD163","F13A1","FCER1G"
)

load("Data/Single_Cell_Data/all_merged_subset_labelled_new.RData")
VlnPlot(all_merged_subset_labelled, 
        features = celltype_markers,
        stack = T,
        alpha = 2,
        sort = F,
        log = F,
        layer = "data",
        fill.by = "ident")+
  NoLegend()+
  labs(x="SCT normalised Log-transformed expression"#,
       #title = "c."
  )+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_line(),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 30),
        strip.text = element_text(size = 9, face = "plain"),
        strip.text.x = element_text(angle = 80, vjust = 0.5,size = 14))

