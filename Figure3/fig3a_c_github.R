library(tidyverse)
library(Seurat)
library(sctransform)

Lv_counts <- read.csv("./orthofinder/Lv6hpf_original.csv")
He_counts <- read.csv("./orthofinder/He6hpf_original.csv")

row.names(He_counts) <- He_counts$Name
He_counts$Name <- NULL

row.names(Lv_counts) <- Lv_counts$Name
Lv_counts$Name <- NULL

set.seed(76)

He_seurat <- CreateSeuratObject(He_counts, project = "he", assay = "RNA",
                                min.cells = 3, min.features = 200, meta.data = NULL)
Lv_seurat <- CreateSeuratObject(Lv_counts, project = "lv", assay = "RNA",
                                min.cells = 3, min.features = 200, meta.data = NULL)

He_seurat$species <- "He"
Lv_seurat$species <- "Lv"

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
ElbowPlot(ortho.transformed.6hpf, ndims = 50) #16
merged_neighbor <- FindNeighbors(object = ortho.transformed.6hpf, dims = 1:10)
merged_clusters <- FindClusters(object = merged_neighbor, resolution = 0.5)
merged_umap_6hpf <- RunUMAP(merged_clusters, dims = 1:10)

merged.markers <- FindAllMarkers(merged_umap_6hpf, 
                                 only.pos = TRUE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)
top_merge <- merged.markers %>% group_by(cluster) %>% top_n(avg_log2FC, n = 40)


merged_umap_6hpf_He <- subset(merged_umap_6hpf, 
                              subset = species == "He" & 
                                seurat_clusters != 10) #removing cluster of putatively dead cells

merged_umap_6hpf_Lv <- subset(merged_umap_6hpf, 
                              subset = species == "Lv")

He_plot <- DimPlot(merged_umap_6hpf_He, reduction = "umap", label = TRUE, pt.size = 0.3) +
  theme(aspect.ratio = 3/4)
Lv_plot <- DimPlot(merged_umap_6hpf_Lv, reduction = "umap", label = TRUE, pt.size = 0.3) +
  theme(aspect.ratio = 3/4)


######################
#Tissue composites
oral_ecto <- list(c("LVA-6953.t1:Sp-Nodal",
                    "LVA-3186.t1:Sp-Gsc",
                    "LVA-28627.t1:Sp-Chordin",
                    "LVA-2999.t1:Sp-Lefty"))
merged_umap_6hpf_Lv <- AddModuleScore(object = merged_umap_6hpf_Lv, 
                                      features = oral_ecto, 
                                      name = "oral_ecto_score")

Lv_oral_ecto <-  FeaturePlot(merged_umap_6hpf_Lv, reduction='umap',
                             features = "oral_ecto_score1",
                             split.by = "species", pt.size = 0.3,
                             order = TRUE)

oral_ecto <- list(c("HER-24103.t1:Sp-Nodal",
                    "HER-18058.t1:Sp-Gsc",
                    "HER-27772.t1:Sp-Chordin",
                    "HER-18254.t1:Sp-Lefty"))
merged_umap_6hpf_He <- AddModuleScore(object = merged_umap_6hpf_He, 
                                      features = oral_ecto, 
                                      name = "oral_ecto_score")

He_oral_ecto <- FeaturePlot(merged_umap_6hpf_He, reduction='umap',
                            features = "oral_ecto_score1",
                            split.by = "species", pt.size = 0.3,
                            order = TRUE)
#endoderm
"LVA-4664.t1:Sp-FoxA",
"LVA-33680.t1:Sp-Blimp1",
"LVA-34456.t1:Sp-Blimp1",
"LVA-34084.t1:Sp-Bra"

endo <- list(c("LVA-4664.t1:Sp-FoxA",
               "LVA-33680.t1:Sp-Blimp1",
               "LVA-34456.t1:Sp-Blimp1",
               "LVA-34084.t1:Sp-Bra"))
merged_umap_6hpf_Lv <- AddModuleScore(object = merged_umap_6hpf_Lv, 
                                      features = endo, 
                                      name = "endo_score")

Lv_endo <- FeaturePlot(merged_umap_6hpf_Lv, reduction='umap',
                       features = "endo_score1",
                       split.by = "species", pt.size = 0.3,
                       order = TRUE)

"HER-17677.t1:Sp-FoxA",
"HER-34700.t1:Sp-Blimp1",
"HER-34701.t1:Sp-Blimp1",
"HER-34716.t1:Sp-Bra",
"HER-1072.t1:Sp-Hnf1-1"

endo <- list(c("HER-17677.t1:Sp-FoxA",
               "HER-34700.t1:Sp-Blimp1",
               "HER-34701.t1:Sp-Blimp1",
               "HER-34716.t1:Sp-Bra",
               "HER-1072.t1:Sp-Hnf1-1"))
merged_umap_6hpf_He <- AddModuleScore(object = merged_umap_6hpf_He, 
                                      features = endo, 
                                      name = "endo_score")

He_endo <- FeaturePlot(merged_umap_6hpf_He, reduction='umap',
                       features = "endo_score1",
                       split.by = "species", pt.size = 0.3,
                       order = TRUE)

#PMC
"LVA-3923.t1:Sp-Delta",
"LVA-5802.t1.2.5f2bf641:Sp-Alx1",
"LVA-5802.t1.3.5f2bf646:Sp-Alx1",
"LVA-6974.t1:Sp-Wnt8",
"LVA-12704.t1:Sp-Tbr"

PMC <- list(c("LVA-3923.t1:Sp-Delta",
              "LVA-5802.t1.2.5f2bf641:Sp-Alx1",
              "LVA-5802.t1.3.5f2bf646:Sp-Alx1",
              "LVA-6974.t1:Sp-Wnt8",
              "LVA-12704.t1:Sp-Tbr"))
merged_umap_6hpf_Lv <- AddModuleScore(object = merged_umap_6hpf_Lv, 
                                      features = PMC, 
                                      name = "pmc_score")

Lv_pmc <- FeaturePlot(merged_umap_6hpf_Lv, reduction='umap',
                      features = "pmc_score1",
                      split.by = "species", pt.size = 0.3,
                      order = TRUE)

"HER-17367.t1:Sp-Delta",
"HER-12794.t1.2.5f2c6949:Sp-Alx1",
"HER-24129.t1:Sp-Wnt8",
"HER-21600.t1:Sp-Tbr"

PMC <- list(c("HER-17367.t1:Sp-Delta",
              "HER-12794.t1.2.5f2c6949:Sp-Alx1",
              "HER-24129.t1:Sp-Wnt8",
              "HER-21600.t1:Sp-Tbr"))
merged_umap_6hpf_He <- AddModuleScore(object = merged_umap_6hpf_He, 
                                      features = PMC, 
                                      name = "pmc_score")

He_pmc <- FeaturePlot(merged_umap_6hpf_He, reduction='umap',
                      features = "pmc_score1",
                      split.by = "species", pt.size = 0.3,
                      order = TRUE)

#APD
"LVA-17097.t1:Sp-FoxQ2-1",
"LVA-17099.t1:Sp-FoxQ2-1",
"LVA-4087.t1:Sp-Six3",
"LVA-4101.t1:Sp-Six3",
"LVA-22122.t1:Sp-Not"

apd <- list(c("LVA-17097.t1:Sp-FoxQ2-1",
              "LVA-17099.t1:Sp-FoxQ2-1",
              "LVA-4087.t1:Sp-Six3",
              "LVA-4101.t1:Sp-Six3",
              "LVA-22122.t1:Sp-Not"))
merged_umap_6hpf_Lv <- AddModuleScore(object = merged_umap_6hpf_Lv, 
                                      features = apd, 
                                      name = "apd_score")

Lv_apd <- FeaturePlot(merged_umap_6hpf_Lv, reduction='umap',
                      features = "apd_score1",
                      split.by = "species", pt.size = 0.3,
                      order = TRUE)

"HER-16550.t1:Sp-FoxQ2-1"
"HER-17866.t1:Sp-Six3"
"HER-9078.t1:Sp-Not"

apd <- list(c("HER-16550.t1:Sp-FoxQ2-1"
              "HER-17866.t1:Sp-Six3"
              "HER-9078.t1:Sp-Not"))
merged_umap_6hpf_He <- AddModuleScore(object = merged_umap_6hpf_He, 
                                      features = apd, 
                                      name = "apd_score")

He_apd <- FeaturePlot(merged_umap_6hpf_He, reduction='umap',
                      features = "apd_score1",
                      split.by = "species", pt.size = 0.3,
                      order = TRUE)
