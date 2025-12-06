library("Seurat")
library("ggplot2")
library("plyr")
library("dplyr")
library("easyGgplot2")
library("ComplexHeatmap")
library("circlize")
library("pcaExplorer")
library("topGO")
library("data.table")
library("tools") 
library("SCENIC")
library("AUCell")
library("SeuratWrappers")
library("monocle3")
library("scales")
library("ggrepel")
library("cowplot")
library("rcartocolor")
library("org.Mm.eg.db")
library("GO.db")
library("readxl")
library("stringr")
library("paletteer")
library("gridExtra")
library("arrow")
library("gghalves")
library("FNN")
library("dplyr")
library("tidyr")
library("purrr")
library("colorspace")
library("tidyverse")

# mLu.DT4 & Xenium
meta.data.list <- list("Cell.Types" = setNames(c("AM", "CD206hi.IMs", "CD206lo.IMs", "Mo-IM", "c.Mo", "nc.Mo", "DC1", "301b+.DC2", "301b−.DC2", "mig.DC", "Neu", "T", "B", "FB", "Endo", "Cyc"), 
                                               c("#AEC7E8FF", "#1F77B4FF", "#DAE3F2", "#F7B6D2FF", "#FB9A4AFF", "#5CB85CFF", "#D95C5CFF", "#A57CBF", "#A9746EFF", "#D798C7FF", "#C3B6D8FF", "#B5BD4BFF", "#53BFD4FF", "#A4CEDDFF", "#C9A18FFF", "#999999FF")))
# mTumor
meta.data.list <- list("Cell.Types" = setNames(c("AMs", "IMs", "Mono", "DC1", "DC2", "Inf.DC2", "Mig.DCs", "DC1-GNs", "GNs", "Eo", "Ba", "MCs", "T8", "T4", "Treg", "Tdn", "Tgd", "NK", "NKP", "B", "pDCs", "ILC2", "B16", "Cycling"), 
                                               c(paletteer_d("ggthemes::Classic_20")[2:1], "#DAE3F2", paletteer_d("ggthemes::Classic_20")[19:20], paletteer_d("ggthemes::Tableau_20")[9:10], "#83AB5B", paletteer_d("ggthemes::Classic_20")[17:18], paletteer_d("ggthemes::Tableau_20")[7:8], paletteer_d("ggthemes::Classic_20")[c(9:10, 13:14)], "#F7E4ED", paletteer_d("ggthemes::Classic_20")[7:8], paletteer_d("ggthemes::Classic_20")[3:4], paletteer_d("ggthemes::Classic_20")[c(11, 5, 15)])))

# mTumor.list
mTumor.HEL1.data <- Read10X(data.dir = "10x_CJ25-40/analysis/CJ_0025_GEX/outs/filtered_feature_bc_matrix")
mTumor.WT1.data <- Read10X(data.dir = "10x_CJ25-40/analysis/CJ_0026_GEX/outs/filtered_feature_bc_matrix")
mTumor.WT2.data <- Read10X(data.dir = "10x_CJ25-40/analysis/CJ_0027_GEX/outs/filtered_feature_bc_matrix")
mTumor.HEL2.data <- Read10X(data.dir = "10x_CJ25-40/analysis/CJ_0028_GEX/outs/filtered_feature_bc_matrix")
mTumor.WT3.data <- Read10X(data.dir = "10x_CJ25-40/analysis/CJ_0033_GEX/outs/filtered_feature_bc_matrix")
mTumor.HEL3.data <- Read10X(data.dir = "10x_CJ25-40/analysis/CJ_0034_GEX/outs/filtered_feature_bc_matrix")
mTumor.HEL.B1.data <- Read10X(data.dir = "10x_CJ25-40/analysis/CJ_0035_GEX/outs/filtered_feature_bc_matrix")
mTumor.muMT1.data <- Read10X(data.dir = "10x_CJ25-40/analysis/CJ_0036_GEX/outs/filtered_feature_bc_matrix")
mTumor.data <- list(mTumor.HEL1.data, mTumor.WT1.data, mTumor.WT2.data, mTumor.HEL2.data, mTumor.WT3.data, mTumor.HEL3.data, mTumor.HEL.B1.data, mTumor.muMT1.data)
mTumor.list.vector <- c("mTumor.HEL1", "mTumor.WT1", "mTumor.WT2", "mTumor.HEL2", "mTumor.WT3", "mTumor.HEL3", "mTumor.HEL.B1", "mTumor.muMT1")
percent.mt.max.y <- 20
nFeature_RNA.min <- c(1500, 1300, 1700, 1500, 600, 1700, 1700, 1950, 1750) # thresholds chosen based on QC plots
nFeature_RNA.max <- c(5000, 5000, 6000, 5000, 2500, 5000, 5000, 5000, 5000)
percent.mt.max <- c(4.8, 5, 5, 6, 4.8, 6, 6.5, 4, 4.8)
mTumor.list <- list()
for (i in 1:length(mTumor.list.vector)){
  mTumor.list[[i]] <- CreateSeuratObject(counts = mTumor.data[[i]], project = mTumor.list.vector[i], min.cells = 3, min.features = 200)
  mTumor.list[[i]]$"orig.ident" <- substr(mTumor.list.vector, 8, nchar(mTumor.list.vector))[i]
  mTumor.list[[i]]$"orig.treatment" <- substr(mTumor.list.vector, 8, nchar(mTumor.list.vector)-1)[i]
  mTumor.list[[i]]$"orig.tissue" <- "mTumor"
  mTumor.list[[i]]$"orig.species" <- "mouse"
  mTumor.list[[i]]$"percent.mt" <- PercentageFeatureSet(mTumor.list[[i]], pattern = "^mt-")
  mTumor.list[[i]] <- subset(mTumor.list[[i]], subset = nFeature_RNA > nFeature_RNA.min[i] & nFeature_RNA < nFeature_RNA.max[i] & percent.mt < percent.mt.max[i])
  mTumor.list[[i]] <- RenameCells(mTumor.list[[i]], add.cell.id = mTumor.list.vector[i])
  mTumor.list[[i]] <- mTumor.list[[i]] %>% SCTransform(vars.to.regress = "percent.mt", method = "glmGamPoi") %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50) %>% FindClusters()
  saveRDS(mTumor.list[[i]], paste0(mTumor.list.vector[i], ".rds"))
}

##mTumor.combined
mTumor.features <- SelectIntegrationFeatures(object.list = mTumor.list, nfeatures = 10000)
mTumor.features.list <- PrepSCTIntegration(object.list = mTumor.list, anchor.features = mTumor.features)
mTumor.anchors <- FindIntegrationAnchors(object.list = mTumor.features.list, normalization.method = "SCT", anchor.features = mTumor.features)
mTumor.combined <- IntegrateData(anchorset = mTumor.anchors, normalization.method = "SCT", new.assay.name = "ITG")
mTumor.combined <- mTumor.combined %>% RunPCA(assay = "ITG") %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50, return.model=TRUE) %>% FindClusters()
# Clusters c(18, 7, 32, 36, 9, 11, 20) are junk cells, expressing all the genes (capturing the environment free mRNA) # remove # redo-cluster
saveRDS(mTumor.combined, paste0("mTumor.combined", ".rds"))

# Pick appropriate resolution for cell type annotation, for increase cluster number
resolution <- c(0.01, 0.02, 0.03, 0.045, 0.05, 0.07, 0.1, 0.12, 0.15, 0.2, 0.22, 0.25, 0.31, 0.32, 0.33, 0.35, 0.39, 0.45, 0.459, 0.47, 0.5, 0.55, 0.6, 0.7, 0.75, 0.81, 0.815, 0.82, 0.86, 0.87, 0.9, 0.92, 0.97, 0.99)
clusters <- vector()
for (i in c(1:length(resolution))){
  mTumor.combined <- FindClusters(mTumor.combined, resolution = resolution[i])
  ITG_snn_res <- paste0("ITG_snn_res.", resolution[i])
  Idents(mTumor.combined) <- ITG_snn_res
  clusters[i] <- length(levels(Idents(mTumor.combined)))
}
names(clusters) <- resolution
clusters
#saveRDS(mTumor.combined, paste0("/Datasets/mTumor.combined.all.res.rds"))
#mTumor.combined <- readRDS(paste0("/Datasets/mTumor.combined.all.res.rds"))
DefaultAssay(mTumor.combined) <- "SCT"
width <- c(seq(from = 1225, by = 125, length.out = 100))
height <- c(seq(from = 270, by = 15, length.out = 100))
for (i in c(1:length(resolution))){
  Idents(mTumor.combined) <- paste0("ITG_snn_res.", resolution[i])
  dimplot <- DimPlot(mTumor.combined, label = TRUE, repel = TRUE, label.size = 6, raster = FALSE)
  mTumor.combined.markers <- FindAllMarkers(mTumor.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  mTumor.combined.markers.top10 <- mTumor.combined.markers %>% group_by(cluster) %>% slice_head(n = 10)
  dotplot <- DotPlot(mTumor.combined, features = unique(mTumor.combined.markers.top10$gene)) + RotatedAxis() + ggtitle(paste0("DEGs_", "ITG_snn_res.", resolution[i])) + ggeasy::easy_center_title()
  dimplot.directory <- paste0("/Datasets/ITG_snn_res/", i+6, ".", "ITG_snn_res.", resolution[i], ".DimPlot.tiff")
  tiff(dimplot.directory, width = 800, height = 640)
  print(dimplot)
  dev.off() 
  dotplot.directory <- paste0("/Datasets/ITG_snn_res/", i+6, ".", "ITG_snn_res.", resolution[i], ".DotPlot.tiff")
  tiff(dotplot.directory, width = width[i], height = height[i])
  print(dotplot)
  dev.off() 
  csv.directory <- paste0("/Datasets/ITG_snn_res/", i+6, ".", "ITG_snn_res.", resolution[i], ".DEGs.csv")
  write.csv(mTumor.combined.markers, csv.directory)
}
#SubCluster for complex cluster
Idents(mTumor.combined) <- "ITG_snn_res.0.33"
DefaultAssay(mTumor.combined) <- "ITG"
mTumor.combined <- FindSubCluster(mTumor.combined, graph.name = "ITG_snn", cluster = 7, resolution = 0.08)

#Cell.Types
# Manual cell-type annotation based on cluster-level marker genes and UMAP neighborhood
# (UMAP coordinate windows here are just a coding convenience to capture coherent regions)
mTumor.combined$Cell.Types <- "Unassigned"
mTumor.combined$Cell.Types[mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < -3 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > -9 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < -5 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > -15] <- "ILC2"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.32 == 19] <- "NKP"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.32 == 8] <- "NK"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.32 == 11] <- "Tgd"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.32 == 3] <- "Tdn"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.32 == 10] <- "T8"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.32 == 5] <- "T4"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.32 == 16] <- "Cycling"
mTumor.combined$Cell.Types[colnames(mTumor.combined) %in% shinyApp.Cell[["T4"]]] <- "T4"
mTumor.combined$Cell.Types[colnames(mTumor.combined) %in% shinyApp.Cell[["Tdn"]]] <- "Tdn"
mTumor.combined$Cell.Types[mTumor.combined$Cell.Types == "T4" & mTumor.combined[["SCT"]]@counts["Foxp3", ] > 0] <- "Treg"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.05 == 3] <- "DC1"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.05 == 4] <- "DC2"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.05 == 5] <- "Mig.DCs"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.05 == 1] <- "IMs"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.32 == 7 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 0 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > -3 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < 0 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > -3] <- "B16"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.33 == 19] <- "Inf.DC2"
mTumor.combined$Cell.Types[mTumor.combined$sub.cluster == "7_1"] <- "Mono"
mTumor.combined$Cell.Types[mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > 13] <- "B" # by marker expression
mTumor.combined$Cell.Types[mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < -13] <- "AMs"
mTumor.combined$Cell.Types[mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < -8 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > -13 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > -3] <- "GNs"
mTumor.combined$Cell.Types[mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < -6 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > -8 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < -3 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > -5] <- "MCs"
mTumor.combined$Cell.Types[mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 13 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > 6 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < 10 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > 0] <- "DC1"
mTumor.combined$Cell.Types[mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 6 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > 5 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < 5 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > 0] <- "DC1-GNs"
mTumor.combined$Cell.Types[mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 5 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > 2 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < 3 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > 0] <- "pDCs"
mTumor.combined$Cell.Types[colnames(mTumor.combined) %in% shinyApp.Cell[["Ba"]]] <- "Ba"
mTumor.combined$Cell.Types[colnames(mTumor.combined) %in% shinyApp.Cell[["Eo"]]] <- "Eo"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.33 %in% c(13, 15, 16, 18)] <- "Cycling"
mTumor.combined$Cell.Types[mTumor.combined$ITG_snn_res.0.7 == 29] <- "Cycling"
mTumor.combined$Cell.Types[mTumor.combined$Cell.Types == "Unassigned" & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 5 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > 1.5 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < 10 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > 7] <- "DC2"
mTumor.combined$Cell.Types[mTumor.combined$Cell.Types == "Unassigned" & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 0 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > -10 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < 20 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > 8] <- "Mig.DCs"
mTumor.combined$Cell.Types[mTumor.combined$Cell.Types == "Unassigned" & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < -2 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > -10 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < 8 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > 0] <- "IMs"
mTumor.combined$Cell.Types[mTumor.combined$Cell.Types == "Unassigned" & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 4 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > -2 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < 7 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > -3] <- "Cycling"
mTumor.combined$Cell.Types[mTumor.combined$Cell.Types == "Unassigned" & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 10 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > 5 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < -2 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > -7] <- "NK"
mTumor.combined$Cell.Types[mTumor.combined$Cell.Types == "Unassigned" & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 10 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > 5 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < -10 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > -20] <- "Tgd"
mTumor.combined$Cell.Types[mTumor.combined$Cell.Types == "Unassigned" & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 4 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > 1 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < -6 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > -11] <- "Tdn"
mTumor.combined$Cell.Types[mTumor.combined$Cell.Types == "Unassigned" & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 3 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > 0.5 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < -2 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > -6] <- "T8"
mTumor.combined$Cell.Types[mTumor.combined$Cell.Types == "Unassigned" & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] < 0.5 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"] > -5 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] < -4 & mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"] > -10] <- "T4"
mTumor.combined$orig.ident <- factor(mTumor.combined$orig.ident, levels = c("WT1", "WT2", "WT3", "muMT1", "HEL1", "HEL2", "HEL3", "HEL.B1"))
mTumor.combined$orig.treatment <- factor(mTumor.combined$orig.treatment, levels = c("WT", "muMT", "HEL", "HEL.B"))
mTumor.combined$Cell.Types <- factor(mTumor.combined$Cell.Types, levels = c("AMs", "IMs", "Mono", "DC1", "DC2", "Inf.DC2", "Mig.DCs", "DC1-GNs", "GNs", "Eo", "Ba", "MCs", "T8", "T4", "Treg", "Tdn", "Tgd", "NK", "NKP", "B", "pDCs", "ILC2", "B16", "Cycling"))
mTumor.combined$Cell.Types[mTumor.combined$seurat_clusters == 8] <- "Mono"
# saveRDS(mTumor.combined, paste0("/Datasets/mTumor.combined.rds"))

####################################################################################################################################################################################

# Figure.1.DimPlot
Figure.1.DimPlot <- DimPlot(mTumor.combined, group.by = "Cell.Types", cols = names(meta.data.list$Cell.Types), raster = FALSE, order = TRUE, pt.size = 1) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range(mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"])*1.05) + ylim(range(mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"])*1.05) + coord_cartesian(expand = FALSE)
Figure.1.DimPlot.directory <- paste0("Figures/Figure.1.DimPlot.pdf")
pdf(Figure.1.DimPlot.directory, width = 18, height = 16)
print(Figure.1.DimPlot)
dev.off() 

# Figure.2.DimPlot
mTumor.combined.Mac <- subset(mTumor.combined, subset = Cell.Types %in% c("IMs", "Mono"))
# saveRDS(mTumor.combined.Mac, paste0("Datasets/SeuratObject/mTumor.combined.Mac.rds"))
Figure.2.DimPlot <- DimPlot(mTumor.combined.Mac, group.by = "Cell.Types", cols = c("#96BADF", "#DAE3F2"), raster = FALSE, order = TRUE, pt.size = 6) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range(mTumor.combined.Mac@reductions$umap@cell.embeddings[, "UMAP_1"])*c(0.95, 1.01)) + ylim(range(mTumor.combined.Mac@reductions$umap@cell.embeddings[, "UMAP_2"])*c(1.01, 0.95)) + coord_cartesian(expand = FALSE)
Figure.2.DimPlot.directory <- paste0("Figures/Figure.2.DimPlot.pdf")
pdf(Figure.2.DimPlot.directory, width = 18, height = 16)
print(Figure.2.DimPlot)
dev.off() 

# Figure.3.FeaturePlot
Figure.3.Feature <- c("Vegfa", "Tgfb1", "Cd274", "Lgals9", "Ceacam1", "Cxcl13", "Tmem119", "Mmp13")
for (i in Figure.3.Feature){
  Figure.3.FeaturePlot <- FeaturePlot(mTumor.combined.Mac, features = i, cols = c("#E5E5E5", "#990F0F"), raster = FALSE, order = TRUE, pt.size = 6) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(range(mTumor.combined.Mac@reductions$umap@cell.embeddings[, "UMAP_1"])*c(0.95, 1.01)) + ylim(range(mTumor.combined.Mac@reductions$umap@cell.embeddings[, "UMAP_2"])*c(1.01, 0.95)) + coord_cartesian(expand = FALSE)
  Figure.3.FeaturePlot.directory <- paste0("Figures/Figure.3.FeaturePlot.", i, ".pdf")
  pdf(Figure.3.FeaturePlot.directory, width = 18, height = 16)
  print(Figure.3.FeaturePlot)
  dev.off() 
}

# Figure.4.VlnPlot
mTumor.combined <- PrepSCTFindMarkers(mTumor.combined)
mTumor.combined.Markers <- FindAllMarkers(mTumor.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(mTumor.combined.Markers, paste0("Datasets/FindAllMarkers/mTumor.combined.Markers.", "DEGs.rds"))
mTumor.combined.Markers.top20 <- mTumor.combined.Markers %>% group_by(cluster) %>% slice_head(n = 20)
features <- c("Ear1", "Ear2", "C1qa", "Csf1r", "Mafb", "Cd14", "Xcr1", "Clec9a", "Mgl2", "Cd209a", "Ifi204", "Ifit2", "Fscn1", "Ccr7", "Cxcr2", "Csf3r", "Adgre1", "Cx3cr1", "Fcer1a", "Cd200r3", "Tpsab1", "Rab27b", "Cd8a", "Nkg7", "Cd4", "Cd5", "Foxp3", "Ctla4", "Tcrg-C1", "Trdc", "Klrk1", "Prf1", "Ms4a1", "Cd19", "Siglech", "Ccr9", "Gata3", "Rora", "Pmel", "Mlana", "Stmn1", "Mki67")
features <- rev(intersect(features, rownames(mTumor.combined[["SCT"]]$data)))
Figure.4.VlnPlot <- VlnPlot(mTumor.combined, group.by = "Cell.Types", stack = TRUE, features = features, fill.by = "ident", cols = names(meta.data.list$Cell.Types), flip = TRUE) + 
  theme(legend.position = "none", axis.title = element_blank(), strip.text.y = element_text(face = "italic", size = 12), axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1.1, size = 18), axis.ticks.x = element_blank())
Figure.4.VlnPlot.directory <- paste0("Figures/Figure.4.VlnPlot.pdf")
pdf(Figure.4.VlnPlot.directory, width = 10, height = 9)
print(Figure.4.VlnPlot)
dev.off()

# Figure.5.DimPlot
mTumor.combined.All.Mac <- mTumor.combined[, mTumor.combined$Cell.Types %in% c("AMs", "IMs", "Mono")]
Figure.5.DimPlot <- DimPlot(mTumor.combined.All.Mac, group.by = "Cell.Types.UMAP", cols = names(meta.data.list$Cell.Types)[1:3], raster = FALSE, order = TRUE, pt.size = 5) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range(mTumor.combined.All.Mac@reductions$umap@cell.embeddings[, "umap_1"])*1.05) + ylim(range(mTumor.combined.All.Mac@reductions$umap@cell.embeddings[, "umap_2"])*1.05) + coord_cartesian(expand = FALSE)
Figure.5.DimPlot.directory <- paste0("Figures/Figure.5.DimPlot.pdf")
pdf(Figure.5.DimPlot.directory, width = 18, height = 16)
print(Figure.5.DimPlot)
dev.off() 

# Figure.6.FeaturePlot
mTumor.combined.Mac.Data <- mTumor.combined.Mac[["SCT"]]@data
mTumor.combined.Mac.Data <- mTumor.combined.Mac.Data[, intersect(colnames(mTumor.combined.Mac.Data), colnames(mTumor.combined.All.Mac))]
mTumor.combined.All.Mac[["SCT"]]$data[, match(colnames(mTumor.combined.Mac.Data), colnames(mTumor.combined.All.Mac))] <- mTumor.combined.Mac.Data
Figure.6.Feature <- c("Mrc1", "Trem2", "Fn1", "Spp1", "Ly6c2", "Pf4", "Plac8", "Metrnl", "Padi4", "Arg1", "Vegfa", "Vcan", "Thbs1", "Tgfb1", "Cd274", "Lgals9", "Ceacam1", "Lyz1", "H2-Ab1", "H2-Aa", "Cxcl13", "Nfe2l2")
for (i in Figure.6.Feature){
  Figure.6.FeaturePlot <- FeaturePlot(mTumor.combined.All.Mac, features = i, cols = c("#E5E5E5", "#990F0F"), raster = FALSE, order = TRUE, pt.size = 5) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(range(mTumor.combined.All.Mac@reductions$umap@cell.embeddings[, "umap_1"])*1.05) + ylim(range(mTumor.combined.All.Mac@reductions$umap@cell.embeddings[, "umap_2"])*1.05) + coord_cartesian(expand = FALSE)
  Figure.6.FeaturePlot.directory <- paste0("Figures/Figure.6.FeaturePlot.", i, ".pdf")
  pdf(Figure.6.FeaturePlot.directory, width = 18, height = 16)
  print(Figure.6.FeaturePlot)
  dev.off() 
}

# Figure.7.DimPlot
Figure.7.DimPlot <- DimPlot(mTumor.combined, group.by = "Cell.Types", cols = names(meta.data.list$Cell.Types), raster = FALSE, order = TRUE, pt.size = 1) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range(mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"])*1.05) + ylim(range(mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"])*1.05) + coord_cartesian(expand = FALSE)
Figure.7.DimPlot.directory <- paste0("Figures/Figure.7.DimPlot.pdf")
pdf(Figure.7.DimPlot.directory, width = 18, height = 16)
print(Figure.7.DimPlot)
dev.off() 

# Figure.8.DimPlot
Figure.8.DimPlot <- DimPlot(mTumor.combined.Mac, group.by = "Cell.Types", cols = names(meta.data.list$Cell.Types)[1:2], raster = FALSE, order = TRUE, pt.size = 5) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range.expand(mTumor.combined.Mac@reductions$umap@cell.embeddings[, "UMAP_1"], 0.02)) + ylim(range.expand(mTumor.combined.Mac@reductions$umap@cell.embeddings[, "UMAP_2"], 0.02)) + coord_cartesian(expand = FALSE)
Figure.8.DimPlot.directory <- paste0("Figures/Figure.8.DimPlot.pdf")
pdf(Figure.8.DimPlot.directory, width = 18, height = 16)
print(Figure.8.DimPlot)
dev.off() 

# Figure.9.FeaturePlot
Figure.9.Feature <- c("Mrc1", "Trem2", "Fn1", "Spp1", "Ly6c2", "Pf4", "Plac8", "Metrnl", "Padi4", "Arg1", "Vegfa", "Vcan", "Thbs1", "Tgfb1", "Cd274", "Lgals9", "Ceacam1", "Lyz1", "H2-Ab1", "H2-Aa", "Cxcl13", "Nfe2l2", "Chil3")
for (i in Figure.9.Feature){
  Figure.9.FeaturePlot <- FeaturePlot(mTumor.combined.Mac, features = i, cols = c("#E5E5E5", "#990F0F"), raster = FALSE, order = TRUE, pt.size = 5) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(range.expand(mTumor.combined.Mac@reductions$umap@cell.embeddings[, "UMAP_1"], 0.02)) + ylim(range.expand(mTumor.combined.Mac@reductions$umap@cell.embeddings[, "UMAP_2"], 0.02)) + coord_cartesian(expand = FALSE)
  Figure.9.FeaturePlot.directory <- paste0("Figures/Figure.9.FeaturePlot.", i, ".pdf")
  pdf(Figure.9.FeaturePlot.directory, width = 18, height = 16)
  print(Figure.9.FeaturePlot)
  dev.off() 
}

# Figure.10.FeaturePlot
Figure.10.Feature <- c("Mrc1", "Trem2", "Fn1", "Spp1", "Ly6c2", "Pf4", "Plac8", "Metrnl", "Padi4", "Arg1", "Vegfa", "Vcan", "Thbs1", "Tgfb1", "Cd274", "Lgals9", "Ceacam1", "Lyz1", "H2-Ab1", "H2-Aa", "Cxcl13", "Nfe2l2", "Chil3")
for (i in Figure.10.Feature){
  Figure.10.FeaturePlot <- FeaturePlot(mTumor.combined, features = i, cols = c("#E5E5E5", "#990F0F"), raster = FALSE, order = TRUE, pt.size = 1) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(range.expand(mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"], 0.05)) + ylim(range.expand(mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"], 0.05)) + coord_cartesian(expand = FALSE)
  Figure.10.FeaturePlot.directory <- paste0("Figures/Figure.10.FeaturePlot.", i, ".pdf")
  pdf(Figure.10.FeaturePlot.directory, width = 18, height = 16)
  print(Figure.10.FeaturePlot)
  dev.off() 
}

# Figure.11.FeaturePlot
Figure.11.Feature <- c("C1qb", "Folr2", "Cd163", "Mmp9", "Lyve1", "Ccr2", "Mmp12", "Cxcl9", "Cxcl10", "C3")
for (i in Figure.11.Feature){
  Figure.11.FeaturePlot <- FeaturePlot(mTumor.combined.Mac, features = i, cols = c("#E5E5E5", "#990F0F"), raster = FALSE, order = TRUE, pt.size = 5) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(range.expand(mTumor.combined.Mac@reductions$umap@cell.embeddings[, "UMAP_1"], 0.02)) + ylim(range.expand(mTumor.combined.Mac@reductions$umap@cell.embeddings[, "UMAP_2"], 0.02)) + coord_cartesian(expand = FALSE)
  Figure.11.FeaturePlot.directory <- paste0("Figures/Figure.11.FeaturePlot.", i, ".pdf")
  pdf(Figure.11.FeaturePlot.directory, width = 18, height = 16)
  print(Figure.11.FeaturePlot)
  dev.off() 
}

# Figure.12.FeaturePlot
Figure.12.Feature <- c("Ccl2", "Ccr2", "Ccl5", "Ccr5", "C1qa", "C1qb", "C1qc", "C3")
for (i in Figure.12.Feature){
  Figure.12.FeaturePlot <- FeaturePlot(mTumor.combined, features = i, cols = c("#E5E5E5", "#990F0F"), raster = FALSE, order = FALSE, pt.size = 2.5) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(range.expand(mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_1"], 0.02)) + ylim(range.expand(mTumor.combined@reductions$umap@cell.embeddings[, "UMAP_2"], 0.02)) + coord_cartesian(expand = FALSE)
  Figure.12.FeaturePlot.directory <- paste0("Figures/Figure.12.FeaturePlot.", i, ".pdf")
  pdf(Figure.12.FeaturePlot.directory, width = 18, height = 16)
  print(Figure.12.FeaturePlot)
  dev.off() 
}

# Figure.13.FeaturePlot
Figure.13.Feature <- c("Mrc1", "Pf4", "Folr2", "Cd163", "Itgax", "Mmp9", "Mmp12")
for (i in Figure.13.Feature){
  Figure.13.FeaturePlot <- FeaturePlot(mLu.combined, features = i, cols = c("#E5E5E5", "#990F0F"), raster = FALSE, order = TRUE, pt.size = 2.5) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
    xlim(range.expand(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_1"], 0.02)) + ylim(range.expand(mLu.combined@reductions$umap@cell.embeddings[, "UMAP_2"], 0.02)) + coord_cartesian(expand = FALSE)
  Figure.13.FeaturePlot.directory <- paste0("Figures/Figure.13.FeaturePlot.", i, ".pdf")
  pdf(Figure.13.FeaturePlot.directory, width = 18, height = 16)
  print(Figure.13.FeaturePlot)
  dev.off() 
}

##Xenium##################################################################################################################################################################################

library("ComplexHeatmap")
library("circlize")
library("viridisLite")
library("grid")

Figure.14.Legend <- Legend(title = "Expression", title_position = "leftcenter-rot", title_gp = gpar(fontsize = 12), at = c(0, 5), labels = c("low", "high"),
                           labels_gp = gpar(fontsize = 14), col_fun = colorRamp2(seq(0, 5, length.out = 100), viridis(100)), legend_height = unit(3, "in"), legend_width = unit(0.6, "in"))
Figure.14.Legend.directory <- paste0("Figures/Figure.14.Legend.pdf")
pdf(Figure.14.Legend.directory, width = 1, height = 4)
print(draw(Figure.14.Legend, x = unit(0.5, "in"), y = unit(2, "in")))
dev.off()

# Figure.15.DotPlot
# Xenium.R # using color #DD0505
# Raw 20240823__201225__Xin_082324/output-XETG00256__0033689__Region_1__20240823__201233/analysis/diffexp/gene_expression_graphclust/differential_expression.csv
Figure.15.DotPlot.Raw.vector <- c("KPAR_1", "KPAR_5", "B16_Big", "B16_Small")
Figure.15.DotPlot.Raw.gene.NO <- rep(c(100, 300), each = 2) # mistake at first with all gene number being 100
for (j in 1:length(Figure.15.DotPlot.Raw.vector)) {
  Figure.15.DotPlot.Raw <- read.csv(paste0("Datasets/Raw/Xenium/diffexp/", Figure.15.DotPlot.Raw.vector[j], "/differential_expression.csv"), row.names = 1)
  row.names(Figure.15.DotPlot.Raw) <- Figure.15.DotPlot.Raw[, 1]
  Figure.15.DotPlot.Raw <- Figure.15.DotPlot.Raw[, -1]
  cluster.NO <- ncol(Figure.15.DotPlot.Raw)/3
  # Figure.15.DotPlot.Data data.frame
  Mean.Counts <- c()
  Log2.fold.change <- c()
  Adjusted.p.value <- c()
  for (i in 1:ncol(Figure.15.DotPlot.Raw)) {
    if ((i+2)%%3 == 0) {Mean.Counts <- c(Mean.Counts, Figure.15.DotPlot.Raw[, i])} else {
      if ((i+2)%%3 == 1) {Log2.fold.change <- c(Log2.fold.change, Figure.15.DotPlot.Raw[, i])} else {
        Adjusted.p.value <- c(Adjusted.p.value, Figure.15.DotPlot.Raw[, i])
      }
    }
  }
  # Cluster.Heatmap
  Heatmap.Data <- Figure.15.DotPlot.Raw[, seq(from = 1, to = ncol(Figure.15.DotPlot.Raw)-2, by = 3)] # Mean.Counts # Decide to use
  colnames(Heatmap.Data) <- paste0("c", str_split_i(colnames(Heatmap.Data), "\\.", 2))
  Heatmap.Data <- scale(t(as.matrix(Heatmap.Data)))
  Cluster.Heatmap <- Heatmap(Heatmap.Data)
  # Figure.15.DotPlot.Data Order
  Figure.15.DotPlot.Data <- data.frame(gene = rep(row.names(Figure.15.DotPlot.Raw), cluster.NO), cluster = paste0("c", rep(1:(ncol(Figure.15.DotPlot.Raw)/3), each = Figure.15.DotPlot.Raw.gene.NO[j])), Mean.Counts = Mean.Counts, Log2.fold.change = Log2.fold.change, Adjusted.p.value = Adjusted.p.value)
  Figure.15.DotPlot.Data$gene <- factor(Figure.15.DotPlot.Data$gene, levels = colnames(Heatmap.Data)[column_order(Cluster.Heatmap)])
  Figure.15.DotPlot.Data$cluster <- factor(Figure.15.DotPlot.Data$cluster, levels = rev(paste0("c", row_order(Cluster.Heatmap))))
  Figure.15.DotPlot.Data$Mean.Counts <- as.vector(t(Heatmap.Data))
  Figure.15.DotPlot.Data$Log2.fold.change[Figure.15.DotPlot.Data$Log2.fold.change < 0] <- 0
  Figure.15.DotPlot.Data$Significant.p.value <- 0
  Figure.15.DotPlot.Data$Significant.p.value[Figure.15.DotPlot.Data$Adjusted.p.value < 0.05 & Figure.15.DotPlot.Data$Log2.fold.change > 0] <- 1
  saveRDS(Figure.15.DotPlot.Data, paste0("Datasets/Figure.Datasets/Figure.15.DotPlot.Data_", Figure.15.DotPlot.Raw.vector[j], ".rds"))
  Figure.15.DotPlot <- ggplot(Figure.15.DotPlot.Data, aes(x = gene, y = cluster, fill = Mean.Counts, size = Log2.fold.change, stroke = Significant.p.value)) + geom_point(shape = 21, color = "black") + scale_fill_gradient(low = "white", high = "#DD0505") + theme_linedraw() + RotatedAxis() +
    theme(axis.title = element_blank(), panel.grid.major = element_line(color = "grey80", size = 0.25, linetype = 2), axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 11, face = "italic"))
  Figure.15.DotPlot.directory <- paste0("Figures/", "Figure.15.DotPlot_", Figure.15.DotPlot.Raw.vector[j], ".pdf")
  pdf(Figure.15.DotPlot.directory, width = ifelse(j %in% 1:2, 23.5, 65.5), height = ifelse(j %in% 1:2, 1.025 + 0.175*cluster.NO, 1.025 + 0.3*cluster.NO))
  print(Figure.15.DotPlot)
  dev.off() 
}

# Figure.16.BarPlot
# Figure.23.BarPlot # for B16F10
mLu.combined.Cell.Types.DEGs <- read.csv("Datasets/FindAllMarkers/Cell.Types.DEGs.Final.csv", row.names = 1)
mLu.combined.IMs.Cell.Types.DEGs <- read.csv("Datasets/FindAllMarkers/IMs.Cell.Types.DEGs.csv", row.names = 1)
Figure.16.BarPlot.Vector <- c("KPAR1", "KPAR5")
Figure.16.BarPlot.Cell.Types <- c("AMs", "CD206hi.IMs", "CD206lo.IMs", "recMacs")
Figure.16.BarPlot.Assignments <- list(c(21, 25, 19, 7), c(14, 19, 17, 6))
Figure.16.BarPlot.Raw <- "Datasets/Raw/Xenium/LassoSelection/"
Figure.16.BarPlot.Raw.Data <- list()
Figure.16.BarPlot.data.frame <- list()
for (i in 1:length(Figure.16.BarPlot.Vector)) {
  Figure.16.BarPlot.Raw.Data[[i]] <- list()
  Figure.16.BarPlot.data.frame[[i]] <- data.frame(row.names = Figure.16.BarPlot.Cell.Types)
  Figure.16.BarPlot.Raw.csv <- list.files(paste0(Figure.16.BarPlot.Raw, Figure.16.BarPlot.Vector[i]))
  for (j in 1:length(Figure.16.BarPlot.Raw.csv)) {
    Figure.16.BarPlot.Raw.Data[[i]][[j]] <- read.csv(paste0(Figure.16.BarPlot.Raw, Figure.16.BarPlot.Vector[i], "/", Figure.16.BarPlot.Raw.csv[j]), skip = 2, header = TRUE)
    Figure.16.BarPlot.Raw.Data[[i]][[j]] <- mapvalues( Figure.16.BarPlot.Raw.Data[[i]][[j]]$Cluster, from = paste0("Cluster ", Figure.16.BarPlot.Assignments[i]), Figure.16.BarPlot.Cell.Types)
    Figure.16.BarPlot.Raw.Data[[i]][[j]] <-  Figure.16.BarPlot.Raw.Data[[i]][[j]][!grepl("Cluster |Unassigned",  Figure.16.BarPlot.Raw.Data[[i]][[j]])]
    Figure.16.BarPlot.table <- table(Figure.16.BarPlot.Raw.Data[[i]][[j]])
    Figure.16.BarPlot.table_full <- setNames(rep(0, length(Figure.16.BarPlot.Cell.Types)), Figure.16.BarPlot.Cell.Types)
    Figure.16.BarPlot.table_full[names(Figure.16.BarPlot.table)] <- Figure.16.BarPlot.table
    Figure.16.BarPlot.table <- Figure.16.BarPlot.table_full
    Figure.16.BarPlot.data.frame[[i]] <- cbind(Figure.16.BarPlot.data.frame[[i]], as.vector(Figure.16.BarPlot.table))
  }
  colnames(Figure.16.BarPlot.data.frame[[i]]) <- paste0("Tumor_",  1:length(Figure.16.BarPlot.Raw.csv))
}
names(Figure.16.BarPlot.data.frame) <- Figure.16.BarPlot.Vector
# saveRDS(Figure.16.BarPlot.data.frame, paste0("Datasets/Figure.Datasets/Figure.16.BarPlot.data.frame.rds"))
for (i in 1:length(Figure.16.BarPlot.Vector)) {
  Figure.16.BarPlot.Data <- as.data.frame(t(Figure.16.BarPlot.data.frame[[i]])/colSums(Figure.16.BarPlot.data.frame[[i]]))
  Figure.16.BarPlot.Data <- melt(Figure.16.BarPlot.Data, varnames = "Cell.Types", value.name = "Freq")
  colnames(Figure.16.BarPlot.Data)[1] <- "Cell.Types"
  Figure.16.BarPlot.Data$orig.ident <- factor(rep(colnames(Figure.16.BarPlot.data.frame[[i]]), 4), levels = rev(as.vector(colnames(Figure.16.BarPlot.data.frame[[i]]))))
  Figure.16.BarPlot.Data$Cell.Types <- factor(Figure.16.BarPlot.Data$Cell.Types, levels = rev(as.vector(levels(Figure.16.BarPlot.Data$Cell.Types))))
  saveRDS(Figure.16.BarPlot.Data, paste0("Datasets/Figure.Datasets/Figure.16.BarPlot.Data_", Figure.16.BarPlot.Vector[i], ".rds"))
  Figure.16.BarPlot <- ggplot2.barplot(Figure.16.BarPlot.Data, xName = "orig.ident", yName= "Freq", width = 0.8, xTickLabelFont = c(22, "plain", color = "grey60"), yTickLabelFont = c(22, "plain", c(rep("black", 4), rep("red", 3))), groupName = "Cell.Types", backgroundColor = "white", removePanelGrid = TRUE, removePanelBorder = TRUE, xShowTitle = FALSE, yShowTitle = FALSE, groupColors = rev(c("#AEC7E8FF", "#1F77B4FF", "#DAE3F2", "#F7B6D2FF"))) + 
    scale_y_continuous(labels = scales::percent) + coord_flip(expand = FALSE) + ylab("Subject") + theme(axis.line = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length = unit(.25, "cm"), axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), legend.position = "none")
  Figure.16.BarPlot.directory <- paste0("Figures/Figure.16.BarPlot_", Figure.16.BarPlot.Vector[i], ".png")
  png(Figure.16.BarPlot.directory, width = 1200, height = 500)
  print(Figure.16.BarPlot)
  dev.off()
}

# Figure.17.DonutChart
reticulate::install_miniconda(force = TRUE)
reticulate::conda_create("r-reticulate", python_version = "3.10")
reticulate::use_condaenv("r-reticulate", required = TRUE)
reticulate::py_install("kaleido", pip = TRUE, envname = "r-reticulate") # to make figures
reticulate::py_install("plotly", pip = TRUE, envname = "r-reticulate") # to make figures
Figure.17.DonutChart.data <- data.frame("KPAR1" = rowSums(Figure.16.BarPlot.data.frame[[1]]), "KPAR5" = rowSums(Figure.16.BarPlot.data.frame[[2]]))
Figure.17.DonutChart.data$Cell.Types <- rownames(Figure.17.DonutChart.data)
reticulate::py_run_string("import sys")
for (i in 1:(length(colnames(Figure.17.DonutChart.data))-1)){
  Figure.17.DonutChart <- Figure.17.DonutChart.data %>% plot_ly(labels = ~Cell.Types, values = ~Figure.17.DonutChart.data[[colnames(Figure.17.DonutChart.data)[i]]], marker = list(colors = c("#AEC7E8FF", "#1F77B4FF", "#DAE3F2", "#F7B6D2FF")), sort = FALSE, direction = "clockwise", automargin = FALSE, textposition = "inside", textfont = list(size = 30))
  Figure.17.DonutChart <- Figure.17.DonutChart %>% add_pie(hole = 0.6) %>% layout(showlegend = F, xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  Figure.17.DonutChart.directory <- paste0("Figures/", "Figure.17.DonutChart.", colnames(Figure.17.DonutChart.data)[i], ".pdf")
  save_image(Figure.17.DonutChart, Figure.17.DonutChart.directory, scale = 5)
}

# Figure.18.ImageDimPlot
home.dir <- "20240823__201225__Xin_082324/"
sample <- c("output-XETG00256__0036733__KPAR_1__20240920__202119", "output-XETG00256__0036733__KPAR_5__20240920__202119")
resolutions <- "gene_expression_graphclust"
sample.name <- c("KPAR_1", "KPAR_2")
Figure.16.BarPlot.Assignments <- list(c(21, 25, 6, 7), c(14, 19, 17, 6))
# mLu.Xenium
for (i in 1:length(sample)) {
  # mLu.Xenium <- LoadXenium(paste0(home.dir, sample[i]), fov = "fov", segmentations = "cell") # long time
  # saveRDS(mLu.Xenium, paste0("Xenium/Datasets/SeuratObject/mLu.Xenium.", sample.name[i], ".UpdatedVersion.rds"))
  clusters <- read.csv(paste0(home.dir, sample[i], "/analysis/clustering/", resolutions, "/clusters.csv"))
  mLu.Xenium$cluster <- clusters$Cluster[match(colnames(mLu.Xenium), clusters$Barcode)]
  mLu.Xenium.Mac <- subset(mLu.Xenium, cluster %in% Figure.16.BarPlot.Assignments[[i]]) # c("#FFFF00", "#FFF700", "#FFEF00", "#FFFF66", "#FFEA00", "#FFD700", "#FFBF00", "#DAA520", "#F5BD1F", "#E3B505", "#FCC200")
  mLu.Xenium.Mac$Cell.Types <- factor(mapvalues(mLu.Xenium.Mac$cluster, from = Figure.16.BarPlot.Assignments[[i]], to = c("AMs", "CD206hi.IMs", "CD206lo.IMs", "recMacs")), levels = c("AMs", "CD206hi.IMs", "CD206lo.IMs", "recMacs"))
  Figure.18.ImageDimPlot <- ImageDimPlot(mLu.Xenium.Mac, fov = "fov", size = 0, molecules = c("Ccl2", "C1qb", "Cd163"), mols.cols = c("#FCC200", "#2350DE", "#990F0F"), mols.size = 0.2, mols.alpha = ifelse(i == 1, 0.2, 0.5), nmols = 20000, dark.background = FALSE, flip_xy = TRUE, axes = FALSE) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")
  Figure.18.ImageDimPlot.directory <- paste0("Figures/Figure.18.ImageDimPlot_", sample.name[i], "_2.png")
  png(Figure.18.ImageDimPlot.directory, width = 600, height = 600)
  grid.newpage()
  grid.arrange(Figure.18.ImageDimPlot, vp = viewport(angle = -90)) # Draw the plot rotated 90° clockwise (i.e. angle = -90)
  dev.off()
}

##B16F10##mMo.DT4##############################################################################################################################################################################

# Figure.19.DimPlot
mMo <- readRDS(paste0("Datasets/BPCells/mMo.Project.RunUMAP.rds"))
mMo.DT4.Raw <- subset(mMo, orig.treatment == "B16")
mMo.DT4 <- subset(mMo.DT4.Raw, Cell.Types.full %in% c("AM", "CD206hi.IMs", "CD206lo.IMs", "Mo-IM", "c.Mo", "nc.Mo", "DC1", "301b+.DC2", "301b−.DC2", "mig.DC", "Neu", "T", "B", "FB", "Endo", "Cyc"))
mMo.DT4$Cell.Types.full <- factor(gsub("Mo-IM", "recMacs", mMo.DT4$Cell.Types.full), levels = c("AM", "CD206hi.IMs", "CD206lo.IMs", "recMacs", "c.Mo", "nc.Mo", "DC1", "301b+.DC2", "301b−.DC2", "mig.DC", "Neu", "T", "B", "FB", "Endo", "Epi", "Cyc"))
mMo.DT4 <- RunUMAP(mMo.DT4, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
# saveRDS(mMo.DT4, "Datasets/SeuratObject/mMo.DT4.rds")
# Figure.19.DimPlot
Figure.19.DimPlot <- DimPlot(mMo.DT4, group.by = "Cell.Types.full", cols = names(meta.data.list$Cell.Types), reduction = "umap.full", raster = FALSE, order = TRUE, pt.size = 1) + theme_void() + theme(plot.title = element_blank()) + 
  xlim(range(mMo.DT4@reductions$umap.full@cell.embeddings[, "UMAPfull_1"])*1.05) + ylim(range(mMo.DT4@reductions$umap.full@cell.embeddings[, "UMAPfull_2"])*1.05) + coord_cartesian(expand = FALSE) +theme(
    legend.text = element_text(size = 40), legend.title = element_text(size = 42), legend.key.size = unit(1, "cm")) + guides(color = guide_legend(override.aes = list(size = 6.5)))
Figure.19.DimPlot.directory <- paste0("Figures/Figure.19.DimPlot.pdf")
pdf(Figure.19.DimPlot.directory, width = 20, height = 16)
print(Figure.19.DimPlot)
dev.off() 

# Figure.20.DotPlot
# Figure.20.DotPlot.DEGs.top3
Idents(mMo.DT4) <- "Cell.Types.full"
mMo.DT4 <- JoinLayers(mMo.DT4)
Figure.20.DotPlot.DEGs <- FindAllMarkers(mMo.DT4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(Figure.20.DotPlot.DEGs, "Datasets/FindAllMarkers/mMo.DT4.DEGs.rds")
Figure.20.DotPlot.DEGs.top3 <- Figure.20.DotPlot.DEGs %>% group_by(cluster) %>% slice_head(n = 3)
Figure.20.DotPlot.DEGs.top3 <- setNames(Figure.20.DotPlot.DEGs.top3$gene, Figure.20.DotPlot.DEGs.top3$cluster)
Figure.20.DotPlot.DEGs.top3 <- Figure.20.DotPlot.DEGs.top3[!duplicated(Figure.20.DotPlot.DEGs.top3)]
Figure.20.DotPlot.DEGs.top3 <- lapply(split(Figure.20.DotPlot.DEGs.top3, names(Figure.20.DotPlot.DEGs.top3)), unname)
Figure.20.DotPlot.DEGs.top3 <- Figure.20.DotPlot.DEGs.top3[order(as.numeric(names(Figure.20.DotPlot.DEGs.top3)))]
saveRDS(Figure.20.DotPlot.DEGs.top3, "Datasets/Figure.Datasets/Figure.20.DotPlot.DEGs.top3.rds")
Figure.20.DotPlot.DEGs.top3 <- list.deduplicate(Figure.20.DotPlot.DEGs.top3)
Figure.20.DotPlot.DEGs.top3 <- Figure.20.DotPlot.DEGs.top3[gsub("Mo-IM", "recMacs", as.vector(meta.data.list$Cell.Types))]
Figure.20.DotPlot <- DotPlot(mMo.DT4, features = Figure.20.DotPlot.DEGs.top3, group.by = "Cell.Types.full", cols = c("white", "#990F0F"), dot.scale = 7.5) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 15), strip.text = element_text(size = 14), legend.position = "none")
Figure.20.DotPlot.directory <- paste0("Figures/Figure.20.DotPlot.DEGs.top3.pdf")
pdf(Figure.20.DotPlot.directory, width = 19, height = 6.8)
print(Figure.20.DotPlot)
dev.off()
# Figure.20.DotPlot.Markers
Figure.20.DotPlot.Markers <- list(
  "AM" = c("Siglecf", "Krt79", "Ear1"),
  "CD206hi.IMs" = c("Folr2", "Cd163", "Pf4"),
  "CD206lo.IMs" = c("C1qb", "Cx3cr1", "H2-Aa"),
  "recMacs" = c("Vcan", "Ccr2", "Ly6c2"),
  "c.Mo" = c("Plac8", "Fn1", "Sell"),
  "nc.Mo" = c("S1pr5", "Ifitm6", "Itgal"),
  "DC1" = c("Itgae", "Xcr1", "Clec9a"),
  "301b+.DC2" = c("Mgl2", "Cd209a", "H2-Oa"),
  "301b−.DC2" = c("Dpp4", "S100a4", "Flt3"),
  "mig.DC" = c("Ccr7", "Ccl22", "Fscn1"),
  "Neu" = c("Cxcr2", "Csf3r", "S100a9"),
  "T" = c("Trbc2", "Cd3g", "Il2rb"),
  "B" = c("Ms4a1", "Igkc", "Cd79a"),
  "FB" = c("Col5a2", "Pdgfra", "Fbn1"),
  "Endo" = c("Cdh5", "Cldn5", "Tek"),
  "Cyc" = c("Top2a", "Stmn1", "Pclaf"))
Figure.20.DotPlot.Markers <- list.deduplicate(Figure.20.DotPlot.Markers)
Figure.20.DotPlot <- DotPlot(mMo.DT4, features = Figure.20.DotPlot.Markers, group.by = "Cell.Types.full", cols = c("white", "#990F0F"), dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 15), strip.text = element_text(size = 14), legend.position = "none")
Figure.20.DotPlot.directory <- paste0("Figures/Figure.20.DotPlot.Markers.pdf")
pdf(Figure.20.DotPlot.directory, width = 19, height = 6.3)
print(Figure.20.DotPlot)
dev.off()

# Figure.21.VlnPlot
dir.create("Figures/Figure.21.VlnPlot")
Figure.21.VlnPlot.Features <- c("C1qb", "Pf4", "Mrc1", "Folr2", "Cd163", "Mmp9", "Lyve1", "Ccr2", "Mmp12", "Fn1", "Ly6c2", "Vcan", "Thbs1", "Spp1", "Vegfa", "Arg1", "Cd274", "Trem2", "Cxcl13", "Cxcl9", "Cxcl10", "Tmem119", "H2-Ab1")
for (i in Figure.21.VlnPlot.Features) {
  Figure.21.VlnPlot <- VlnPlot(mMo.DT4, group.by = "Cell.Types.full", idents = c("CD206hi.IMs", "CD206lo.IMs", "recMacs"), features = i, cols = names(meta.data.list$Cell.Types)[2:4]) +
    theme(axis.title.x = element_blank(), axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 10), legend.position = "") + ggtitle(NULL)
  Figure.21.VlnPlot.directory <- paste0("Figures/Figure.21.VlnPlot/Figure.21.VlnPlot_", i, ".pdf")
  pdf(Figure.21.VlnPlot.directory, width = 4, height = 3.5)
  print(Figure.21.VlnPlot)
  dev.off()
}

# Figure.22.FeaturePlot
DefaultAssay(mMo.DT4) <- "RNA"
Figure.22.FeaturePlot <- FeaturePlot(mMo.DT4, features = "Ccl2", reduction = "umap.full", raster = FALSE, order = TRUE, pt.size = 0.75) + scale_colour_gradient(low = "#E5E5E5", high = "#990F0F") + 
  theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range.expand(mMo.DT4@reductions$umap.full@cell.embeddings[, "UMAPfull_1"], expand.ratio = 0.05)) + ylim(range.expand(mMo.DT4@reductions$umap.full@cell.embeddings[, "UMAPfull_2"], expand.ratio = 0.05)) + coord_cartesian(expand = FALSE)
Figure.22.FeaturePlot.directory <- paste0("Figures/Figure.22.FeaturePlot.png")
png(Figure.22.FeaturePlot.directory, width = 1100, height = 1000)
print(Figure.22.FeaturePlot)
dev.off()

# Figure.24.DonutChart
reticulate::install_miniconda(force = TRUE)
reticulate::conda_create("r-reticulate", python_version = "3.10")
reticulate::use_condaenv("r-reticulate", required = TRUE)
reticulate::py_install("kaleido", pip = TRUE, envname = "r-reticulate") # to make figures
reticulate::py_install("plotly", pip = TRUE, envname = "r-reticulate") # to make figures
Figure.24.DonutChart.data <- data.frame("B16_Big" = rowSums(Figure.23.BarPlot.data.frame[[1]]), "B16_Small" = rowSums(Figure.23.BarPlot.data.frame[[2]]))
Figure.24.DonutChart.data$Cell.Types <- rownames(Figure.24.DonutChart.data)
reticulate::py_run_string("import sys")
for (i in 1:(length(colnames(Figure.24.DonutChart.data))-1)){
  Figure.24.DonutChart <- Figure.24.DonutChart.data %>% plot_ly(labels = ~Cell.Types, values = ~Figure.24.DonutChart.data[[colnames(Figure.24.DonutChart.data)[i]]], marker = list(colors = c("#AEC7E8FF", "#1F77B4FF", "#DAE3F2", "#F7B6D2FF")), sort = FALSE, direction = "clockwise", automargin = FALSE, textposition = "inside", textfont = list(size = 30))
  Figure.24.DonutChart <- Figure.24.DonutChart %>% add_pie(hole = 0.6) %>% layout(showlegend = F, xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  Figure.24.DonutChart.directory <- paste0("Figures/", "Figure.24.DonutChart.", colnames(Figure.24.DonutChart.data)[i], ".pdf")
  save_image(Figure.24.DonutChart, Figure.24.DonutChart.directory, scale = 5)
}

# Figure.25.ImageDimPlot
# To fix arrow error form v22.0.0 "Error: Invalid: Invalid number of indices: 0"; might get fixed in future updates
# See transcripts.csv.gz.R script to make transcripts.csv.gz
home.dir <- "250307_Xenium/20250225__211706__Xin_022525/"
sample <- c("output-XETG00256__0043844__BG16_Big__20250225__211728", "output-XETG00256__0043844__BG16_Small__20250225__211728")
sample.name <- c("B16F10_1", "B16F10_2")
Cell.Types <- c("AMs", "CD206hi.IMs", "CD206lo.IMs", "recMacs")
Cell.Types.Assignments <- list(list(15, c(33, 18), c(25, 28, 24, 46), c(36, 40, 41, 24, 26, 23, 20)), list(c(15, 27, 43), c(21, 11, 23), c(35, 25, 26, 20), c(31, 37, 39, 18, 41, 14)))
resolutions <- "gene_expression_graphclust"
for (i in 1:length(sample)) {
  # mLu.Xenium <- LoadXenium(paste0(home.dir, sample[i]), fov = "fov", segmentations = "cell") # long time
  # saveRDS(mLu.Xenium, paste0("Xenium/Datasets/SeuratObject/mLu.Xenium.", sample.name[i], ".UpdatedVersion.rds"))
  clusters <- read.csv(paste0(home.dir, sample[i], "/analysis/clustering/", resolutions, "/clusters.csv"))
  mLu.Xenium$cluster <- clusters$Cluster[match(colnames(mLu.Xenium), clusters$Barcode)]
  mLu.Xenium$Cell.Types <- factor("Others", levels = c("AMs", "CD206hi.IMs", "CD206lo.IMs", "recMacs"))
  for (j in 1:4) {
    mLu.Xenium$Cell.Types[mLu.Xenium$cluster %in% Cell.Types.Assignments[[i]][[j]]] <- Cell.Types[j]
  }
  mLu.Xenium.Mac <- subset(mLu.Xenium, Cell.Types %in% Cell.Types)
  Figure.25.ImageDimPlot <- ImageDimPlot(mLu.Xenium.Mac, fov = "fov", size = 0, molecules = c("Ccl2", "C1qb", "Cd163"), mols.cols = c("#FFD700", "#2E6AD3", "#990F0F"), mols.size = 0.2, mols.alpha = 0.2, nmols = 20000, dark.background = FALSE, flip_xy = TRUE, axes = FALSE) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")
  Figure.25.ImageDimPlot.directory <- paste0("Figures/Figure.25.ImageDimPlot_", sample.name[i], ".png")
  png(Figure.25.ImageDimPlot.directory, width = 600, height = 600)
  grid.newpage()
  grid.arrange(Figure.25.ImageDimPlot, vp = viewport(angle = -90)) # Draw the plot rotated 90° clockwise (i.e. angle = -90)
  dev.off()
}

# Figure.26.ImageDimPlot
# Figure.27.ImageDimPlot
# Figure.28.ImageDimPlot
ImageDimPlot.directory.vector <- c("Figures/Figure.26.ImageDimPlot_Stromal/", 
                                   "Figures/Figure.27.ImageDimPlot_Myeloid/", 
                                   "Figures/Figure.28.ImageDimPlot_Markers/")
sample <- c("20240823__201225__Xin_082324/output-XETG00256__0036733__KPAR_1__20240920__202119", "20240823__201225__Xin_082324/output-XETG00256__0036733__KPAR_5__20240920__202119", 
            "250307_Xenium/20250225__211706__Xin_022525/output-XETG00256__0043844__BG16_Big__20250225__211728", "250307_Xenium/20250225__211706__Xin_022525/output-XETG00256__0043844__BG16_Small__20250225__211728")
resolutions <- "gene_expression_graphclust"
sample.name <- c("KPAR_1", "KPAR_2", "B16F10_1", "B16F10_2")
Cell.Types.Assignments <- list(
  KPAR_1 = list(AMs = 21, CD206hi.IMs = 25, CD206lo.IMs = 6, recMacs = 7,
                Tumor = c(11, 17, 22, 24), # Tumor # Based on expression pattern
                Epi = c(8), # Epi # Epcam
                Endo = c(18)), # Endo # Based on expression pattern
  KPAR_2 = list(AMs = 14, CD206hi.IMs = 19, CD206lo.IMs = 17, recMacs = 6,
                Tumor = c(22, 5), # Tumor # Based on expression pattern
                Epi = c(8), # Epi # Epcam
                Endo = c(12)), # Endo # Based on expression pattern
  B16F10_1 = list(AMs = 15, CD206hi.IMs = c(33, 18), CD206lo.IMs = c(25, 28, 24, 46), recMacs = c(36, 40, 41, 24, 26, 23, 20),
                  Tumor = c(10, 8, 31, 16, 22), # Tumor # Pmel
                  Epi = c(47, 42, 6), # Epi # Epcam
                  Endo = c(19)), # Endo # Pecam1, Kdr, Esm1, Eng, Tek
  B16F10_2 = list(AMs = c(15, 27, 43), CD206hi.IMs = c(21, 11, 23), CD206lo.IMs = c(35, 25, 26, 20), recMacs = c(31, 37, 39, 18, 41, 14),
                  Tumor = c(12, 19, 30), # Tumor # Pmel
                  Epi = c(5), # Epi # Epcam
                  Endo = c(13))) # Endo # Pecam1, Kdr, Esm1, Eng, Tek
# mLu.Xenium <- list()
for (i in 1:length(sample)) {
  mLu.Xenium[[i]] <- readRDS(paste0("Xenium/Datasets/SeuratObject/mLu.Xenium.", sample.name[i], ".UpdatedVersion.rds"))
  clusters <- read.csv(paste0(sample[i], "/analysis/clustering/", resolutions, "/clusters.csv"))
  mLu.Xenium[[i]]$cluster <- clusters$Cluster[match(colnames(mLu.Xenium[[i]]), clusters$Barcode)]
  mLu.Xenium[[i]]$Cell.Types <- as.vector("Others")
  for (j in 1:length(Cell.Types.Assignments[[i]])) {
    mLu.Xenium[[i]]$Cell.Types[mLu.Xenium[[i]]$cluster %in% Cell.Types.Assignments[[i]][[j]]] <- names(Cell.Types.Assignments[[i]])[j]
  }
  mLu.Xenium[[i]]$Cell.Types <- factor(mLu.Xenium[[i]]$Cell.Types, levels = c(names(Cell.Types.Assignments[[1]]), "Others"))
  saveRDS(mLu.Xenium[[i]], paste0("Xenium/Datasets/SeuratObject/mLu.Xenium.", sample.name[i], ".UpdatedVersion.rds"))
  Figure.26.ImageDimPlot <- ImageDimPlot(mLu.Xenium[[i]], cells = colnames(mLu.Xenium[[i]])[mLu.Xenium[[i]]$Cell.Types %in% names(Cell.Types.Assignments[[1]])[5:7]], group.by = "Cell.Types", cols = c("grey70", names(meta.data.list$Cell.Types)[c(8, 10)]), fov = "fov", size = 1, dark.background = FALSE, flip_xy = TRUE, axes = FALSE) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")
  Figure.27.ImageDimPlot <- ImageDimPlot(mLu.Xenium[[i]], cells = colnames(mLu.Xenium[[i]])[mLu.Xenium[[i]]$Cell.Types %in% names(Cell.Types.Assignments[[1]])[1:4]], group.by = "Cell.Types", cols = names(meta.data.list$Cell.Types)[c(15, 2, 1, 4)], fov = "fov", size = 1, dark.background = FALSE, flip_xy = TRUE, axes = FALSE) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")
  Figure.28.ImageDimPlot <- ImageDimPlot(mLu.Xenium[[i]], fov = "fov", size = 0, molecules = c("Cxcl9", "Cxcl10", "Cxcl13", "Ccl2"), mols.cols = c(gsub("#DF5C24FF", "#F8B323", paletteer_d("ggthemes::few_Dark")[c(3, 9, 4, 2)])), mols.size = 0.2, mols.alpha = 0.2, nmols = 20000, dark.background = FALSE, flip_xy = TRUE, axes = FALSE) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + theme(legend.position = "none")
  ImageDimPlot.list <- list(Figure.26.ImageDimPlot, Figure.27.ImageDimPlot, Figure.28.ImageDimPlot)
  for (k in 1:3) {
    dir.create(ImageDimPlot.directory.vector[k], recursive = TRUE, showWarnings = FALSE)
    ImageDimPlot.directory <- paste0(ImageDimPlot.directory.vector[k], sample.name[i], "_", c("Stromal", "Myeloid", "Markers")[k], ".png")
    png(ImageDimPlot.directory, width = 600, height = 600)
    grid.newpage()
    grid.arrange(ImageDimPlot.list[[k]], vp = viewport(angle = -90)) # Draw the plot rotated 90° clockwise (i.e. angle = -90)
    dev.off()
  }
}

# Figure.29.VlnPlot
# distances in pixel units; plotted as log1p(distance)
dir.create("Figures/Figure.29.VlnPlot/", recursive = TRUE, showWarnings = FALSE)
for (i in 1:length(sample)) {
  cell_coords <- mLu.Xenium[[i]]@images$fov$centroids@coords
  coords_AMs <- cell_coords[mLu.Xenium[[i]]$Cell.Types == "AMs", ]
  coords_CD206hi.IMs <- cell_coords[mLu.Xenium[[i]]$Cell.Types == "CD206hi.IMs", ]
  coords_CD206lo.IMs <- cell_coords[mLu.Xenium[[i]]$Cell.Types == "CD206lo.IMs", ]
  coords_recMacs <- cell_coords[mLu.Xenium[[i]]$Cell.Types == "recMacs", ]
  coords_Tumor <- cell_coords[mLu.Xenium[[i]]$Cell.Types == "Tumor", ]
  coords_BV <- cell_coords[mLu.Xenium[[i]]$Cell.Types %in% c("Epi", "Endo"), ]
  # Compute distances from each Myeloid Cells to Stromal Cells
  AMs_to_Tumor <- get.knnx(data = coords_AMs, query = coords_Tumor, k = 1)$nn.dist
  CD206hi.IMs_to_Tumor <- get.knnx(data = coords_CD206hi.IMs, query = coords_Tumor, k = 1)$nn.dist
  CD206lo.IMs_to_Tumor <- get.knnx(data = coords_CD206lo.IMs, query = coords_Tumor, k = 1)$nn.dist
  recMacs_to_Tumor <- get.knnx(data = coords_recMacs, query = coords_Tumor, k = 1)$nn.dist
  AMs_to_BV <- get.knnx(data = coords_AMs, query = coords_BV, k = 1)$nn.dist
  CD206hi.IMs_to_BV <- get.knnx(data = coords_CD206hi.IMs, query = coords_BV, k = 1)$nn.dist
  CD206lo.IMs_to_BV <- get.knnx(data = coords_CD206lo.IMs, query = coords_BV, k = 1)$nn.dist
  recMacs_to_BV <- get.knnx(data = coords_recMacs, query = coords_BV, k = 1)$nn.dist
  # ---- Put everything into one data frame ----
  df <- tibble(
    value  = c(AMs_to_Tumor, CD206hi.IMs_to_Tumor, CD206lo.IMs_to_Tumor, recMacs_to_Tumor,
               AMs_to_BV, CD206hi.IMs_to_BV, CD206lo.IMs_to_BV, recMacs_to_BV),
    group  = factor(rep(c("T1","T2","T3","T4","B1","B2","B3","B4"),
                        times = c(length(AMs_to_Tumor), length(CD206hi.IMs_to_Tumor), length(CD206lo.IMs_to_Tumor), length(recMacs_to_Tumor),
                                  length(AMs_to_BV), length(CD206hi.IMs_to_BV), length(CD206lo.IMs_to_BV), length(recMacs_to_BV))),
                    levels = c("T1","T2","T3","T4","B1","B2","B3","B4")),
    type   = factor(rep(c("Tumor","Tumor","Tumor","Tumor",
                          "Bronchovascular bundle","Bronchovascular bundle","Bronchovascular bundle","Bronchovascular bundle"),
                        times = c(length(AMs_to_Tumor), length(CD206hi.IMs_to_Tumor), length(CD206lo.IMs_to_Tumor), length(recMacs_to_Tumor),
                                  length(AMs_to_BV), length(CD206hi.IMs_to_BV), length(CD206lo.IMs_to_BV), length(recMacs_to_BV))))
  )
  # ---- Log-transform ----
  df <- df %>% mutate(value_log = log1p(value))
  # ---- Plot with log scale ----
  Figure.29.VlnPlot <- ggplot(df, aes(x = group, y = value_log, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "black") +
    geom_boxplot(width = 0.15, outlier.size = 0.5, color = "black") +
    facet_grid(. ~ type, scales = "free_x", space = "free_x") +
    scale_fill_brewer(palette = "Set2") +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank(), strip.text = element_text(color = "black", size = 14), axis.text.x = element_text(angle = 45, hjust = 1, size = 14), legend.position = "none"
    ) + labs(x = "", y = "log1p(Distance)"
    ) + scale_x_discrete(labels = c("T1" = "AMs", "T2" = "CD206hi.IMs", "T3" = "CD206lo.IMs", "T4" = "recMacs", "B1" = "AMs", "B2" = "CD206hi.IMs", "B3" = "CD206lo.IMs", "B4" = "recMacs"
    )) + scale_fill_manual(values = setNames(rep(names(meta.data.list$Cell.Types)[c(15, 2, 1, 4)], 2), paste0(rep(c("T", "B"), each = 4), c(1:4, 1:4))))
  Figure.29.VlnPlot.directory <- paste0("Figures/Figure.29.VlnPlot/Figure.29.VlnPlot_", sample.name[i], ".pdf")
  pdf(Figure.29.VlnPlot.directory, width = 6, height = 5)
  print(Figure.29.VlnPlot)
  dev.off()
}

