# Script to complete initial clustering of the scHuBMAP data into stromal, immune, and epithelial cells

parent_directory <- "/hubmap_single_cell/"

# Import libraries
.libPaths("./libraries/R_LIBS_4p1p2/")
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(stringr)
library(DoubletFinder)
library(harmony)
#library(celda)
`%notin%` <- Negate(`%in%`)
set.seed(1)

# load helper functions
source(paste0(parent_directory, "scRNA/scripts/single_cell_rna_helper_functions.R"))

# load the metadata (only needed for step 2)
metadata <- read.table(paste0(parent_directory, "/hubmap_metadata_atac_and_rna_final.csv"), header = TRUE, sep = ",", stringsAsFactors=FALSE)

# set steps to run
execute_steps <- c(1,2,3,4,5,6,7)

# Define variables
sample_name <- "all_samples"

# folder containing project after mergining individual projects post QC
individual_qc_and_dublet_plot_location <- paste0(parent_directory, "scRNA/initial_processing/")

# folder to save results
analysis_parent_folder <- paste0(parent_directory, "scRNA/all_cells/decontx_norm_scale_harmony/")

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 1) Add decontX counts to the seurat object and set as default assay
if (1 %in% execute_steps){
	# load initial project
	setwd(individual_qc_and_dublet_plot_location)
	colon <- readRDS("full_intestine_proj_seurat.rds")
	colon[["percent.ribo"]] <- PercentageFeatureSet(colon, pattern = "^RP[SL]")

	# Add decontx counts
	full_counts <- readRDS(paste0(parent_directory, "scRNA/initial_processing/fulldecontxobj.rds"))
	decontXcounts <- CreateAssayObject(counts = GetAssayData(object = full_counts, slot = "counts")) # create a new assay to store counts
	colon[["decontXcounts"]] <- decontXcounts

	# Set decontx as default assay
	DefaultAssay(colon) <- "decontXcounts"

	# set directory to store new results
	setwd(analysis_parent_folder)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 2) Add metadata
if (2 %in% execute_steps){
	metadata <- metadata[,colnames(metadata)[2:dim(metadata)[2]]]
	colnames(metadata) <- c("Sample", colnames(metadata)[2:dim(metadata)[2]])
	metadata <- metadata[metadata$Sample != "",]

	meta_data_types <- colnames(metadata)
	for (i in 2:length(meta_data_types)){
		identities <- colon[['orig.ident']]
		for (j in 1:length(metadata$Sample)){
			identities[identities==metadata$Sample[j]] <- metadata[j,meta_data_types[i]]
		}
		colon <- AddMetaData(colon, identities$orig.ident, col.name = meta_data_types[i])
	}
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 3) Normalize and scale data
if (3 %in% execute_steps){
	colon <- seurat_standard_normalize_and_scale(colon, TRUE, 1.0)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 4) Plot UMAPs and run harmony
if (4 %in% execute_steps){
	# Plot by clustering, sample, and disease state
	pdf(paste0("./UMAPclustering" , ".pdf"), onefile=F)
	print(DimPlot(colon, reduction = "umap", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_samples.pdf"), width = 12, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))#, cols = (ArchRPalettes$stallion))
	dev.off()

	pdf(paste0("./UMAP_location.pdf"), width = 12, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "Location", cols = c("#208A42", "#272E6A", "#D51F26", "#FEE500", "#F8AD16", "#F47D2B", "#89288F", "#438496")))#, cols = (ArchRPalettes$stallion))
	dev.off()

	colon <- RunHarmony(colon, "Donor") # will use PCA
	colon <- RunUMAP(colon, dims = 1:20, reduction = "harmony", reduction.name = "umapharmony")
	pdf(paste0("./UMAPharmony_clustering.pdf"), width = 12)
	print(DimPlot(colon, reduction = "umapharmony", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAPharmony_samples.pdf"), width = 12)
	print(DimPlot(colon, reduction = "umapharmony", group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))#, cols = (ArchRPalettes$stallion))
	dev.off()

	pdf(paste0("./UMAPharmony_location.pdf"), width = 12, onefile=F)
	print(DimPlot(colon, reduction = "umapharmony", group.by = "Location", cols = c("#208A42", "#272E6A", "#D51F26", "#FEE500", "#F8AD16", "#F47D2B", "#89288F", "#438496")))#, cols = (ArchRPalettes$stallion))
	dev.off()

	colon <- FindNeighbors(colon, reduction = "harmony", dims = 1:20)
	set.seed(1)
	colon <- FindClusters(colon, resolution = 0.5)

	pdf(paste0("./UMAPharmony_newclustering.pdf"), width = 12)
	plot = DimPlot(colon, reduction = "umapharmony", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE))
	print(LabelClusters(plot = plot, id = "ident"))
	dev.off()

	saveRDS(colon, "clustered_full_colon_proj_seurat.rds")
}


##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 5) Plot Markers
if (5 %in% execute_steps){
	reductions_to_plot <- c("umapharmony")
	for (reduction in reductions_to_plot){
		seurat_feature_plot(colon, sample_name, reduction, "B_cells", c("PAX5", "MS4A1", "CD19", "IGLL5", "VPREB3"))
		seurat_feature_plot(colon, sample_name, reduction, "GC_B_cells", c("SERPINA9", "HRK", "HTR3A", "TCL6", "CD180", "FCRLA"))
		seurat_feature_plot(colon, sample_name, reduction, "Plasma_B_cells", c("SSR4", "IGLL5", "IGLL1", "AMPD1"))
		seurat_feature_plot(colon, sample_name, reduction, "Mast_cells", c("TPSAB1", "HDC", "CTSG", "CMA1", "KRT1", "IL1RAPL1", "GATA2"))
		seurat_feature_plot(colon, sample_name, reduction, "Immature_Goblet", c("KLK1","ITLN1","WFDC2","CLCA1","LRRC26","RETNLB","SPINK4","AGR2"))
		seurat_feature_plot(colon, sample_name, reduction, "Secretory", c("MUC6", "MUC5B"))
		seurat_feature_plot(colon, sample_name, reduction, "Liver cells", c("ALB", "BHMT", "TAT", "CPS1"))
		seurat_feature_plot(colon, sample_name, reduction, "EnteroendocrineMarkers", c("CCK", "GIP", "PYY", "INSL5", "SST", "MLN", "TPH1", "PROX1", "SMOC2", "RPX3", "RPX6", "ONECUT3", "GHRL", "GCG"))
		seurat_feature_plot(colon, sample_name, reduction, "CyclingTA", c("TICRR","CDC25C","SAPCD2","MT3","LINC00669","SERPINA6"))
		seurat_feature_plot(colon, sample_name, reduction, "Secretory", c("MUC6", "MUC5B"))
		seurat_feature_plot(colon, sample_name, reduction, "Goblet", c("MUC2", "TFF1", "FCGBP","FFAR4","SYTL2","LGALS9B","BCAS1"))
		seurat_feature_plot(colon, sample_name, reduction, "Stem", c("SMOC2", "RGMB", "LGR5", "ASCL2", "SOX9", "CD34"))
		seurat_feature_plot(colon, sample_name, reduction, "Enteroendocrine", c("CRYBA2","SCGN","FEV","CHGA","GCG","SCG5","PCSK1N","PYY","NEUROD1","MS4A8","DDC"))
		seurat_feature_plot(colon, sample_name, reduction, "Tuft", c("GNG13","SH2D7","SH2D6","TRPM5","AZGP1","KRT18","BMX","PSTPIP2","LRMP","PTGS1","IL17RB","HCK","PLCG2","ANXA13"))
		seurat_feature_plot(colon, sample_name, reduction, "Fibroblasts", c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "FAP", "CBLN2", "SPOCK1", "ACSS3"))
		seurat_feature_plot(colon, sample_name, reduction, "Fibroblasts_cell_subtypes", c("RSPO3", "CCL11", "WNT5B", "BMP4", "CHI3L1", "ACTA2", "WNT2B"))
		seurat_feature_plot(colon, sample_name, reduction, "Myofibroblasts", c("SYT10", "SOSTDC1", "DES", "TAGLN"))
		seurat_feature_plot(colon, sample_name, reduction, "stromal_other", c("FAM110D", "INHBB", "NPR1", "NOVA2", "GPIHBP1", "SOX17", "VWF", "PLVAP", "CDH5", "S100B"))
		seurat_feature_plot(colon, sample_name, reduction, "CD69pos_Mast", c("CMA1", "IL1RAPL1", "CD69"))
		seurat_feature_plot(colon, sample_name, reduction, "NK", c("KLRF1", "SH2D1B", "SH2D1B"))
		seurat_feature_plot(colon, sample_name, reduction, "Monocytes_macrophages", c("CD14", "CLEC9A", "FCGR1A", "LILRB2", "CD209", "CD1E", "FOLR2","FABP3","PLA2G2D"))
		seurat_feature_plot(colon, sample_name, reduction, "Enterocytes", c("CA1", "RAB6B"))
		seurat_feature_plot(colon, sample_name, reduction, "Best4posEnterocytes", c("BEST4", "CA7","OTOP2","OTOP3", "MYOM1","MT1G","MT1H"))
		seurat_feature_plot(colon, sample_name, reduction, "General_Epithelial", c("EPCAM", "KRT8","KRT18"))
		seurat_feature_plot(colon, sample_name, reduction, "T_cells", c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "TBX21", "IL7R", "CD4", "CD2"))
		seurat_feature_plot(colon, sample_name, reduction, "Tregs", c("BATF","TNFRSF4", "FOXP3","CTLA4","LAIR2"))
	}
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 6) ID Cluster Markers
if (6 %in% execute_steps){
	downsampled.obj <- colon[, sample(colnames(colon), size = 40000, replace=F)]

	# Find cluster specific markers
	colon.markers <- FindAllMarkers(downsampled.obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = 250)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 7) Initial cluster identification
if (7 %in% execute_steps){
	colon <- readRDS("clustered_full_colon_proj_seurat.rds")
	# Rough initial cluster identification
	new.cluster.ids <- c("Epithelial", #0
		"Epithelial",#
		"Epithelial",# 
		"Epithelial",#
		"Stromal",#
		"Immune",#5
		"Epithelial",# 
		"Epithelial",# 
		"Stromal",# 
		"Immune",#
		"Stromal",#10
		"Epithelial",#
		"Epithelial",#
		"Epithelial",#
		"Epithelial",#
		"Immune",#15
		"Stromal",#
		"Epithelial",#
		"Stromal",#
		"Immune",# 
		"Stromal",#20 
		"Epithelial",#
		"Epithelial",#
		"Stromal",#
		"Epithelial",#
		"Stromal",#25
		"Epithelial",# 
		"Epithelial",#
		"Epithelial", #
		"Epithelial", #
		"Epithelial",#30
		"Epithelial", #
		"Immune", #32
		"Stromal", #
		"Epithelial", #
		"Immune",#35
		"Epithelial",#
		"Stromal",#
		"Immune", #
		"Epithelial", #
		"Epithelial",#40
		"Epithelial",#
		"Epithelial",#
		"Epithelial",#
		"Stromal",#
		"Stromal",#45
		"Immune",#46
		"Epithelial"#
		)

	identities <- as.character(colon@meta.data$seurat_clusters)
	for (i in 0:length(new.cluster.ids)){
		identities[identities==as.character(i)] <- new.cluster.ids[i+1]
	}
	colon <- AddMetaData(colon, identities, col.name = "CellTypeInitial")
	pdf(paste0("./umapharmony_cell_type_initial.pdf"), width = 8)
	AugmentPlot(DimPlot(colon, reduction = "umapharmony", group.by = "CellTypeInitial", cols = paletteDiscrete(values = unique(colon@meta.data$CellTypeInitial), set = "stallion", reverse = FALSE))+theme_ArchR())#, cols = (ArchRPalettes$stallion))
	dev.off()

	pdf(paste0("./umapharmony_cell_type_initial_v3.pdf"), width = 8)
	plot = DimPlot(colon, reduction = "umapharmony", group.by = "CellTypeInitial", cols = paletteDiscrete(values = unique(colon@meta.data$CellTypeInitial), set = "stallion", reverse = FALSE))
	plot = LabelClusters(plot = plot, id = "CellTypeInitial")
	AugmentPlot(plot)
	dev.off()

	# Save ids of immune, epithelial, and stromal cells
	immune_clusters <- which(new.cluster.ids == "Immune")-1
	stromal_clusters <- which(new.cluster.ids == "Stromal")-1
	epithelial_clusters <- which(new.cluster.ids == "Epithelial")-1

	# Subset to make immune, stromal, and epithelial projects
	colon <- DietSeurat(colon)
	colon_immune <- DietSeurat(subset(colon, subset = seurat_clusters %in% immune_clusters))
	saveRDS(colon_immune, file = "colon_immune_all_samples_initial.rds")

	colon_stromal <- DietSeurat(subset(colon, subset = seurat_clusters %in% stromal_clusters))
	saveRDS(colon_stromal, file = "colon_stromal_all_samples_initial.rds")

	colon_epitehlial <- DietSeurat(subset(colon, subset = seurat_clusters %in% epithelial_clusters))
	saveRDS(colon_epitehlial, file = "colon_epithelial_all_samples_initial.rds")
}

