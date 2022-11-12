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
analysis_parent_folder <- paste0(parent_folder, "scRNA/epithelial/enteroendocrine/decontx_sctransform_cca/") # folder for this analysis
setwd(analysis_parent_folder)

# load helper functions
source(paste0(parent_folder, "scRNA_reproduce/scripts/single_cell_rna_helper_functions.R"))

# Load and merge data
colon <- readRDS(parent_folder, "scRNA/epithelial/colon/decontx_sctransform_cca/colon_clustered.rds"))
ileum <- readRDS(parent_folder, "scRNA/epithelial/ileum/decontx_sctransform_cca/ileum.rds"))
jejunum <- readRDS(parent_folder, "scRNA/epithelial/jejunum/decontx_sctransform_cca/jejunum_clustered.rds"))
duodenum <- readRDS(parent_folder, "scRNA/epithelial/duodenum/decontx_sctransform_cca/duodenum_clustered.rds"))
intestine <- merge(colon, y = c(ileum, jejunum, duodenum), project = "full_intestine_project")

# Split datasets and run sctransform
DefaultAssay(intestine) <- "decontXcounts"
intestine <- DietSeurat(
    intestine,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    features = NULL,
    assays = c("RNA", "decontXcounts"),
    dimreducs = NULL,
    graphs = NULL
    )
intestine <- DietSeurat(intestine)
rna.list <- SplitObject(object = intestine , split.by = "orig.ident")
for (i in names(rna.list)){
  rna.list[[i]] <- SCTransform(rna.list[[i]], assay = "decontXcounts", method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo"), verbose = TRUE)
}

#integrate everything with the following references
# B001-A-001 46
# B001-A-301 1
# B005-A-002 10
# B006-A-501 71
# B008-A-301 45
features <- SelectIntegrationFeatures(object.list = rna.list, nfeatures = 3000)
rna.list <- PrepSCTIntegration(object.list = rna.list, anchor.features = features)
rna.list <- lapply(X = rna.list, FUN = RunPCA, features = features)
rna.anchors <- FindIntegrationAnchors(object.list = rna.list, normalization.method = "SCT",
    anchor.features = features, dims = 1:30, reference = c(1,10,45,46,71))
rna.combined.sct <- IntegrateData(anchorset = rna.anchors, normalization.method = "SCT", dims = 1:30)

# Now subset to only the enteroendocrine cells
rna.combined.sct.new <- subset(rna.combined.sct, subset = CellType %in% c("Enteroendocrine"))
rna.combined.sct.new <- RunPCA(rna.combined.sct.new, verbose = TRUE)
rna.combined.sct.new <- RunUMAP(rna.combined.sct.new, reduction = "pca", dims = 1:30)
plotUMAP(rna.combined.sct.new, reduction = "umap", save_name = "rna_combined_sct_post_filter1")
set.seed(1)
rna.combined.sct.new <- FindNeighbors(rna.combined.sct.new, reduction = "pca", dims = 1:30)
rna.combined.sct.new <- FindClusters(rna.combined.sct.new, resolution = 3.0)
plotUMAP(rna.combined.sct.new, reduction = "umap", save_name = "rna_combined_sct_ref1_new_clust_post_filter1")

# Plot markers
DefaultAssay(rna.combined.sct.new) <- "decontXcounts"
reduction <- c("umap")
seurat_feature_plot(rna.combined.sct.new, "enteroendocrine_integrated", reduction, "Enteroendocrine", c("CHGA", "TPH1", "CPE", "SST", "CCK", "ONECUT3", "GIP", "MLN", "GHRL", "SCT", "PROX1", "SMOC2", "GCG", "PYY", "INSL5"))

# Find markers
enteroendocrine.markers <- FindAllMarkers(rna.combined.sct.new, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = 250)
head(enteroendocrine.markers[enteroendocrine.markers$cluster == 6,],20)

# Set the cell types in the project
DefaultAssay(rna.combined.sct.new) <- "integrated"
new.cluster.ids <- c("Enterochromaffin", #0
        "L Cells", #1
        "EnteroendocrineUn", #2-
        "EnteroendocrineUn", #3-
        "EnteroendocrineUn 1", #4
        "I Cells", #5
        "Mo Cells", #6
        "EnteroendocrineUn", #7-
        "EnteroendocrineUn", #8-
        "Enterochromaffin", #9
        "Enterochromaffin", #10
        "EnteroendocrineUn", #11
        "EnteroendocrineUn 1", #12
        "D Cells", #13
        "K Cells", #14
        "EnteroendocrineUn", #15
        "Enterochromaffin", #16
        "NEUROG3high", #17
        "S Cells", #18
        "I Cells", #19
        "EnteroendocrineUn")#20

identities <- data.frame(rna.combined.sct.new[['seurat_clusters']])
identities$seurat_clusters <- as.character(rna.combined.sct.new[['seurat_clusters']]$seurat_clusters)
for (i in 0:(length(new.cluster.ids)-1)){
  identities[identities$seurat_clusters==as.character(i),] <- new.cluster.ids[i+1]
}

# plot
rna.combined.sct.new <- AddMetaData(rna.combined.sct.new, identities$seurat_clusters, col.name = "CellType")
pdf(paste0("./enteroendocrine_UMAP_cell_type_id.pdf"), width = 6, onefile = FALSE)
DimPlot(rna.combined.sct.new, group.by = "CellType", reduction = "umap", cols = paletteDiscrete(values = unique(rna.combined.sct.new@meta.data$CellType), set = "kelly", reverse = FALSE)) + theme_ArchR()
dev.off()

# Save
saveRDS(rna.combined.sct.new, "clustered_annotated_enteroendocrine_cells.rds")

