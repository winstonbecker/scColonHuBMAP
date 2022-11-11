# Script to analyze stromal scRNA cells from ENCODE colon data
# WRB 2022

# Import libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ArchR)
library(viridis)
library(DoubletFinder)
library(Rcpp)
library(harmony)
library(future)
library(Matrix)
library(wesanderson)
library(MAST)
library(ggpubr)
library(harmony)
library(SingleR)
library(celldex)
set.seed(1)

# Define variables
parent_directory <- "/hubmap_single_cell/" #folder for the project
sample_name <- "colon" # used as label for saving plots
analysis_parent_folder <- paste0(parent_folder, "scRNA_reproduce/epithelial/colon/decontx_sctransform_cca/") # folder for this analysis
setwd(analysis_parent_folder)

# load helper functions
source(paste0(parent_folder, "scRNA_reproduce/scripts/single_cell_rna_helper_functions.R"))

# Load previously defined seurat object that contains only the cells classified as immune cells
ileum <- readRDS(paste0(parent_folder, "scRNA_reproduce/all_cells/decontx_norm_scale_harmony/colon_epithelial_samples_initial.rds"))

# Markers from Regev Lab paper
RNA_markers_epithelial <- read.csv(file = paste0(parent_folder, "RegevLabMarkers/S2_Regev_Cell_epithelial.csv"))
RNA_markers_immune <- read.csv(file = paste0(parent_folder, "RegevLabMarkers/S2_Regev_Cell_immune.csv"))
RNA_markers_stromal <- read.csv(file = paste0(parent_folder, "RegevLabMarkers/S2_Regev_Cell_stromal.csv"))
RNA_markers <- rbind(RNA_markers_epithelial, RNA_markers_immune, RNA_markers_stromal)

set.seed(1)

###############################################################################################################################
###############################################################################################################################
# SCT and integrate approach
colon[["percent.ribo"]] <- PercentageFeatureSet(colon, pattern = "^RP[SL]")

rna.list <- SplitObject(object = colon , split.by = "orig.ident")
for (i in names(rna.list)){
  rna.list[[i]] <- SCTransform(rna.list[[i]], assay = "decontXcounts", method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo"), verbose = TRUE)
}

features <- SelectIntegrationFeatures(object.list = rna.list, nfeatures = 3000)
rna.list <- PrepSCTIntegration(object.list = rna.list, anchor.features = features)
rna.list <- lapply(X = rna.list, FUN = RunPCA, features = features)
rna.anchors <- FindIntegrationAnchors(object.list = rna.list, normalization.method = "SCT",
    anchor.features = features, dims = 1:30, reference = c(1,20)) #reduction = "rpca", k.anchor = 20, 
rna.combined.sct <- IntegrateData(anchorset = rna.anchors, normalization.method = "SCT", dims = 1:30)

saveRDS(rna.combined.sct, file = "colon_initial_integration.rds")
rna.combined.sct <- readRDS("colon_initial_integration.rds")

set.seed(1)
rna.combined.sct <- RunPCA(rna.combined.sct, verbose = TRUE)
rna.combined.sct <- RunUMAP(rna.combined.sct, reduction = "pca", dims = 1:30)
plotUMAP(rna.combined.sct, reduction = "umap", save_name = "rna_combined_sct_ref1")
set.seed(1)
rna.combined.sct <- FindNeighbors(rna.combined.sct, reduction = "pca", dims = 1:30)
rna.combined.sct <- FindClusters(rna.combined.sct, resolution = 0.5)
plotUMAP(rna.combined.sct, reduction = "umap", save_name = "rna_combined_sct_ref1_new_clust")

# rna.combined.sct <- find_anchors_label_transfer(rna.combined.sct, 
#   seRNA_path = "/oak/stanford/groups/wjg/wbecker/other/scRNA/hubmap_analysis/epithelial_cells/colon/colon_epithelial_cells_metadata_cleanup.rds",
#   label = "hubmap_colon", 
#   reduction = "umap", pred_name = "CellType")

DefaultAssay(rna.combined.sct) <- "decontXcounts"
reductions_to_plot <- c("umap")
for (reduction in reductions_to_plot){
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "CyclingTA", c("TICRR","CDC25C","SAPCD2","MT3","LINC00669","SERPINA6"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "ImmatureEnterocytes", c("SLC26A2","CA1", "KIAA1239","GLRA4","CELA3A","CELA3B","C8G"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "Mcells", c("SOX8","NTRK2","CCL23","CCL20","CTSE","SLC2A6","TNFAIP2","POLD1","GJB3","KCNE2","TNFRSF11A","AKR1C2","SPINK5","ATP11B","NOXO1","SDC4"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "Tuft", c("GNG13","SH2D7","SH2D6","TRPM5","AZGP1","KRT18","BMX","PSTPIP2","LRMP","PTGS1","IL17RB","HCK","PLCG2","ANXA13"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "Enterocytes", c("CA1", "RAB6B"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "Best4posEnterocytes", c("BEST4", "CA7","OTOP2","OTOP3", "MYOM1","MT1G","MT1H"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "General_Epithelial", c("EPCAM", "KRT8","KRT18"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "Immature_Goblet", c("KLK1","ITLN1","WFDC2","CLCA1","LRRC26","RETNLB","SPINK4","AGR2"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "Goblet", c("MUC2", "TFF1", "FCGBP","FFAR4","SYTL2","LGALS9B","BCAS1"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "Stem", c("SMOC2", "RGMB", "LGR5", "ASCL2", "SOX9", "CD34", "BMI1"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "Enteroendocrine", c("CRYBA2","SCGN","FEV","CHGA","GCG","SCG5","PCSK1N","PYY","NEUROD1","MS4A8","DDC"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "Paneth", c("LYZ", "DEFA5", "DEFA6", "REG3A", "PRSS2", "PLA2G2A", "ITLN2", "GUCA2A", "TFF3", "IFI6"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "M cells", c("GP2", "KRT7", "NTRK2","CCL23","CCL20","POLD1","GJB3","KCNE2","RNF207","TNFRSF11A","AKR1C2","SPINK5","NOXO1"))
  seurat_feature_plot(rna.combined.sct, sample_name, reduction, "Other cells", c("CD3", "CD4", "CD8", "CD14", "MS4A4", "PLVAP", "VWF")) 
}
DefaultAssay(rna.combined.sct) <- "decontXcounts"

rna.combined.sct.markers <- FindAllMarkers(rna.combined.sct, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = 250)
RNA_markers_spec <- RNA_markers[RNA_markers$either,]
RNA_markers_spec[RNA_markers_spec$gene %in% rna.combined.sct.markers[rna.combined.sct.markers$cluster == 1,]$gene,"ident"]

DefaultAssay(rna.combined.sct) <- "integrated"

pdf(paste0("./UMAP_pMT_R1.pdf"), width = 12, onefile=F)
print(FeaturePlot(rna.combined.sct, reduction = "umap", features = "percent.mt"))#, cols = (ArchRPalettes$stallion))
dev.off()

pdf(paste0("./UMAP_pRIBO_R1.pdf"), width = 12, onefile=F)
print(FeaturePlot(rna.combined.sct, reduction = "umap", features = "percent.ribo"))#, cols = (ArchRPalettes$stallion))
dev.off()


# remove bad clusters and redo
bad_clusters <- c(1,10,11)
set.seed(1)
# Subset
rna.combined.sct.new <- subset(rna.combined.sct, subset = seurat_clusters %ni% bad_clusters)
rna.combined.sct.new <- RunPCA(rna.combined.sct.new, verbose = TRUE)
rna.combined.sct.new <- RunUMAP(rna.combined.sct.new, reduction = "pca", dims = 1:30)
plotUMAP(rna.combined.sct.new, reduction = "umap", save_name = "rna_combined_sct_post_filter1")
set.seed(1)
rna.combined.sct.new <- FindNeighbors(rna.combined.sct.new, reduction = "pca", dims = 1:30)
rna.combined.sct.new <- FindClusters(rna.combined.sct.new, resolution = 2.0)
plotUMAP(rna.combined.sct.new, reduction = "umap", save_name = "rna_combined_sct_ref1_new_clust_post_filter1")

DefaultAssay(rna.combined.sct.new) <- "decontXcounts"
reductions_to_plot <- c("umap")
for (reduction in reductions_to_plot){
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "CyclingTA", c("TICRR","CDC25C","SAPCD2","MT3","LINC00669","SERPINA6"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "ImmatureEnterocytes", c("SLC26A2","CA1", "KIAA1239","GLRA4","CELA3A","CELA3B","C8G"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "Mcells", c("SOX8","NTRK2","CCL23","CCL20","CTSE","SLC2A6","TNFAIP2","POLD1","GJB3","KCNE2","TNFRSF11A","AKR1C2","SPINK5","ATP11B","NOXO1","SDC4"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "Tuft", c("GNG13","SH2D7","SH2D6","TRPM5","AZGP1","KRT18","BMX","PSTPIP2","LRMP","PTGS1","IL17RB","HCK","PLCG2","ANXA13"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "Enterocytes", c("CA1", "RAB6B"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "Best4posEnterocytes", c("BEST4", "CA7","OTOP2","OTOP3", "MYOM1","MT1G","MT1H"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "General_Epithelial", c("EPCAM", "KRT8","KRT18"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "Immature_Goblet", c("KLK1","ITLN1","WFDC2","CLCA1","LRRC26","RETNLB","SPINK4","AGR2"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "Goblet", c("MUC2", "TFF1", "FCGBP","FFAR4","SYTL2","LGALS9B","BCAS1"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "Stem", c("SMOC2", "RGMB", "LGR5", "ASCL2", "SOX9", "CD34", "BMI1"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "Enteroendocrine", c("CRYBA2","SCGN","FEV","CHGA","GCG","SCG5","PCSK1N","PYY","NEUROD1","MS4A8","DDC"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "Paneth", c("LYZ", "DEFA5", "DEFA6", "REG3A", "PRSS2", "PLA2G2A", "ITLN2", "GUCA2A", "TFF3", "IFI6"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "M cells", c("GP2", "KRT7", "NTRK2","CCL23","CCL20","POLD1","GJB3","KCNE2","RNF207","TNFRSF11A","AKR1C2","SPINK5","NOXO1"))
  seurat_feature_plot(rna.combined.sct.new, sample_name, reduction, "Other cells", c("CD3", "CD4", "CD8", "CD14", "MS4A4", "PLVAP", "VWF")) 
}

rna.combined.sct.markers <- FindAllMarkers(rna.combined.sct.new, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = 250)
RNA_markers_spec <- RNA_markers[RNA_markers$either,]
RNA_markers_spec[RNA_markers_spec$gene %in% rna.combined.sct.markers[rna.combined.sct.markers$cluster == 1,]$gene,"ident"]

DefaultAssay(rna.combined.sct.new) <- "integrated"

# rna.combined.sct.new <- find_anchors_label_transfer(rna.combined.sct.new, 
#   seRNA_path = "/oak/stanford/groups/wjg/wbecker/other/hubmap_single_cell/scRNA/epithelial/colon/counts_sctransform_cca/colon_clustered.rds",
#   label = "hubmap_colon", 
#   reduction = "umap", pred_name = "CellType")

new.cluster.ids <- c("TA2",
        "TA1",
        "TA1",
        "TA2",
        "TA1",
        "CyclingTA 1",
        "Stem",
        "Immature Goblet",
        "Immature Enterocytes",
        "Immature Enterocytes",
        "Enterocytes",
        "Immature Enterocytes",
        "TA1",
        "CyclingTA 1",
        "CyclingTA 1",
        "Best4+ Enterocytes",
        "Goblet",
        "Enterocytes",
        "Immature Goblet",
        "Enterocytes",
        "Immature Goblet",
        "Immature Goblet",
        "Best4+ Enterocytes",
        "CyclingTA 2",
        "Enteroendocrine",
        "Enterocytes",
        "CyclingTA 1",
        "Tuft")#27

# Set the cell types in the project
identities <- data.frame(rna.combined.sct.new[['seurat_clusters']])
identities$seurat_clusters <- as.character(rna.combined.sct.new[['seurat_clusters']]$seurat_clusters)
for (i in 0:(length(new.cluster.ids)-1)){
  identities[identities$seurat_clusters==as.character(i),] <- new.cluster.ids[i+1]
}

pal <- c("#371377")
names(pal) <- "Best4+ Enterocytes"
pal["CyclingTA 1"] <- "grey35"
pal["CyclingTA 2"] <- "#cccccc"
pal["Enterocytes"] <- "#9E0142"
pal["Enteroendocrine"] <- "#FF0080"
pal["Goblet"] <- "#DC494C"
pal["Immature Enterocytes"] <- "#7700FF"
pal["Immature Goblet"] <- "#FAD510"
pal["Stem"] <- "#02401B"
pal["TA1"] <- "#079290"
pal["TA2"] <- "#A2A475"
pal["Tuft"] <- "#0AD7D3"

# plot
rna.combined.sct.new <- AddMetaData(rna.combined.sct.new, identities$seurat_clusters, col.name = "CellType")
pdf(paste0("./UMAP_cell_type_id.pdf"), width = 6, onefile = FALSE)
DimPlot(rna.combined.sct.new, group.by = "CellType", reduction = "umap", cols = pal) + theme_ArchR()
dev.off()

saveRDS(rna.combined.sct.new, file = "colon_clustered.rds")

