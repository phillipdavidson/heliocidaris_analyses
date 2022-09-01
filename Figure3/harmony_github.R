library(tidyverse)
library(Seurat)
library(sctransform)
library(harmony)
set.seed(76)

He_counts <- read.csv("~/Desktop/scRNAseq/He6hpf_integrated.csv")
row.names(He_counts) <- He_counts$Name
He_counts$Name <- NULL

He_seurat <- CreateSeuratObject(He_counts, project = "he", assay = "RNA",
                                min.cells = 3, min.features = 200, meta.data = NULL)

Lv_counts <- read.csv("~/Desktop/scRNAseq/Lv6hpf_integrated.csv")
row.names(Lv_counts) <- Lv_counts$Name
Lv_counts$Name <- NULL

Lv_seurat <- CreateSeuratObject(Lv_counts, project = "lv", assay = "RNA",
                                min.cells = 3, min.features = 200, meta.data = NULL)
He_seurat$group <- "He6hpf"
Lv_seurat$group <- "Lv6hpf"

He_seurat$species <- "He"
Lv_seurat$species <- "Lv"

He_seurat$timepoint <- "6hpf"
Lv_seurat$timepoint <- "6hpf"

ortho.merged.6hpf <- merge(Lv_seurat, y = c(He_seurat), 
                           add.cell.ids = c("Lv6hpf","He6hpf"), 
                           project = "ortho_merged")
ortho.merged.6hpf <- subset(ortho.merged.6hpf, 
                            subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 2500 & 
                              nCount_RNA < 7500)

ortho.transformed.6hpf <- SCTransform(ortho.merged.6hpf, 
                                      variable.features.n = 3000, 
                                      verbose = FALSE)

ortho.transformed.6hpf <- RunPCA(object = ortho.transformed.6hpf, npcs = 50, 
                                 features = VariableFeatures(object = ortho.transformed.6hpf))
merged_neighbor <- FindNeighbors(object = ortho.transformed.6hpf, dims = 1:21)
merged_clusters <- FindClusters(object = merged_neighbor, resolution = 0.5)
merged_umap_6hpf <- RunUMAP(merged_clusters, dims = 1:21)

merged.markers <- FindAllMarkers(merged_umap_6hpf, 
                                 only.pos = TRUE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)

merged.list <- SplitObject(merged_umap_6hpf, split.by = "species")
merged.list <- lapply(X = merged.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = merged.list, nfeatures = 3000)
merged.list <- PrepSCTIntegration(object.list = merged.list, 
                                  anchor.features = features)
merged.anchors <- FindIntegrationAnchors(object.list = merged.list, normalization.method = "SCT",
                                         anchor.features = features)
merged.combined.sct <- IntegrateData(anchorset = merged.anchors, normalization.method = "SCT")
merged.combined.sct <- RunPCA(merged.combined.sct, verbose = FALSE, npcs = 50)
Ortho6hpf <- RunUMAP(merged.combined.sct, reduction = "pca", dims = 1:21)

combined_harmony <- merged.combined.sct %>%
  RunHarmony("species", assay.use="SCT", dims = 1:50)

# plot variance explained by PC:
ElbowPlot(combined_harmony, ndims = 50)

combined_harmony <- RunUMAP(combined_harmony, reduction='harmony', dims = 1:21)
combined_harmony <- FindNeighbors(combined_harmony, reduction='harmony')

harmony_sp_plot <- 
  DimPlot(combined_harmony, group.by='species', reduction='umap', label = TRUE, 
          repel = TRUE, cols = c("darkorange2","seagreen4"), pt.size = 0.3) +
  theme(aspect.ratio = 3/4)

harmony_split_plot <- 
  DimPlot(combined_harmony, group.by='seurat_clusters', split.by = "species",
          reduction='umap', label = TRUE, 
          repel = TRUE, pt.size = 0.3) +
  theme(aspect.ratio = 3/4)
