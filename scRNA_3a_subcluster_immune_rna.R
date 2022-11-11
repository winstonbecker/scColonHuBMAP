# Script to analyze immune scRNA cells from HuBMAP colon data
# WRB 2022

# 1) Initial dim reduciton and clustering
# 2) Iteratively remove bad clusters and redo analysis
# 3) Cluster identification

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
parent_directory <- "/hubmap_single_cell/"
sample_name <- "all_samples" # used as label for saving plots
analysis_parent_folder <- paste0(parent_folder, "scRNA_reproduce/immune/decontx_norm_scale_harmony/") # folder for this analysis
setwd(analysis_parent_folder)

# load helper functions
source(paste0(parent_folder, "scRNA_reproduce/scripts/single_cell_rna_helper_functions.R"))

# Set steps to run
execute_steps <- c(1,2,3)

# Load previously defined seurat object that contains only the cells classified as immune cells
colon <- readRDS(paste0(parent_folder, "scRNA_reproduce/all_cells/decontx_norm_scale_harmony/colon_immune_all_samples_initial.rds"))

# Markers from Regev Lab paper
RNA_markers_epithelial <- read.csv(file = paste0(parent_folder, "RegevLabMarkers/S2_Regev_Cell_epithelial.csv"))
RNA_markers_immune <- read.csv(file = paste0(parent_folder, "RegevLabMarkers/S2_Regev_Cell_immune.csv"))
RNA_markers_stromal <- read.csv(file = paste0(parent_folder, "RegevLabMarkers/S2_Regev_Cell_stromal.csv"))
RNA_markers <- rbind(RNA_markers_epithelial, RNA_markers_immune, RNA_markers_stromal)

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 1) Initial dim reduciton and clustering
if (1 %in% execute_steps){
	colon <- seurat_standard_normalize_and_scale(colon, cluster = FALSE, ndims <- 20)
	plotUMAP(colon, reduction = "umap", save_name = "pre_filter")

	pdf(paste0("./UMAP_location.pdf"), width = 12, onefile=F)
	print(DimPlot(colon, reduction = "umap", group.by = "Location", cols = c("#208A42", "#272E6A", "#D51F26", "#FEE500", "#F8AD16", "#F47D2B", "#89288F", "#438496")))#, cols = (ArchRPalettes$stallion))
	dev.off()

	colon <- run_harmony(colon)
	set.seed(1)
	colon <- FindClusters(colon, resolution = 0.5)
	plotUMAP(colon, reduction = "umapharmony", save_name = "pre_filter")

	pdf(paste0("./UMAPharmony_location.pdf"), width = 12, onefile=F)
	print(DimPlot(colon, reduction = "umapharmony", group.by = "Location", cols = c("#208A42", "#272E6A", "#D51F26", "#FEE500", "#F8AD16", "#F47D2B", "#89288F", "#438496")))#, cols = (ArchRPalettes$stallion))
	dev.off()

	colon.markers <- FindAllMarkers(colon, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = 250)

	# Look at overlapping specific markers with previously defined marker set and print out overlapping cases
	RNA_markers_spec <- RNA_markers[RNA_markers$either,]
	for (val in unique(colon.markers$cluster)){
	    clusterSpecficMarkers <- RNA_markers_spec[RNA_markers_spec$gene %in% colon.markers[colon.markers$cluster == val,]$gene,]
	    print(val)
	    print(clusterSpecficMarkers$ident)
	}

	colon <- run_singleR_and_plot(colon, types_to_use  = c("DC","Epithelial_cells","B_cell","Neutrophils","T_cells","Monocyte",
		"Endothelial_cells","Neurons","Macrophage","NK_cell",
		"BM","Platelets","Fibroblasts","Astrocyte","Myelocyte","Pre-B_cell_CD34-","Pro-B_cell_CD34+","Pro-Myelocyte"))

	# colon <- find_anchors_label_transfer(colon, 
	# 	seRNA_path = "htan_immune_object.rds",
	# 	label = "htanimmune", 
	# 	reduction = "umapharmony")

	# We will move possible doublet clusters
	# Save ids of immune, epithelial, and stromal cells
	stromal_clusters <- c(8)
	epithelial_clusters <- c(7,16)

	# Subset to make immune, stromal, and epithelial projects
	colon_stromal <- DietSeurat(subset(colon, subset = seurat_clusters %in% stromal_clusters))
	colon_epitehlial <- DietSeurat(subset(colon, subset = seurat_clusters %in% epithelial_clusters))

	saveRDS(colon_epitehlial, file = "immune_inital_move_to_epithelail.rds")
	saveRDS(colon_stromal, file = "immune_inital_move_to_stromal.rds")
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 2) Iteratively remove bad clusters and redo analysis
if (2 %in% execute_steps){
	# Save ids of cells to remove
	bad_clusters <- c(8,7,16)

	# Subset seurat object
	colon_new <- DietSeurat(subset(colon, subset = seurat_clusters %ni% bad_clusters))
	sample_name <- paste0(sample_name, "doublet_removed")
	#colon_new <- sc_transform_umap_cluster(colon_new)
	colon_new <- seurat_standard_normalize_and_scale(colon_new, cluster = FALSE, ndims <- 20)
	plotUMAP(colon_new, reduction = "umap", save_name = "post_filter_1")

	colon <- colon_new
	set.seed(1)
	colon <- run_harmony(colon, resolution = 1.8)
	set.seed(1)
	colon <- FindClusters(colon, resolution = 1.2)
	plotUMAP(colon, reduction = "umapharmony", save_name = "post_filter_1")

	reductions_to_plot <- c("umapharmony")
	for (reduction in reductions_to_plot){
		seurat_feature_plot(colon, sample_name, reduction, "B_cells", c("PAX5", "MS4A1", "CD19", "IGLL5", "VPREB3"))
		seurat_feature_plot(colon, sample_name, reduction, "GC_B_cells", c("SERPINA9", "HRK", "HTR3A", "TCL6", "CD180", "FCRLA"))
		seurat_feature_plot(colon, sample_name, reduction, "Plasma_B_cells", c("SSR4", "IGLL5", "IGLL1", "AMPD1"))
		seurat_feature_plot(colon, sample_name, reduction, "Mast_cells", c("TPSAB1", "HDC", "CTSG", "CMA1", "KRT1", "IL1RAPL1", "GATA2"))
		seurat_feature_plot(colon, sample_name, reduction, "CD69pos_Mast", c("CMA1", "IL1RAPL1", "CD69"))
		seurat_feature_plot(colon, sample_name, reduction, "NK", c("KLRF1", "SH2D1B", "SH2D1B", "NCAM1", "FCGR3A")) # NCAM1 and FCGR3A from Blish lab paper
		seurat_feature_plot(colon, sample_name, reduction, "Monocytes_macrophages", c("CD14", "CLEC9A", "FCGR1A", "LILRB2", "CD209", "CD1E", "FOLR2","FABP3","PLA2G2D"))
		seurat_feature_plot(colon, sample_name, reduction, "T_cells", c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "TBX21", "IL7R", "CD4", "CD2"))
		seurat_feature_plot(colon, sample_name, reduction, "Tregs", c("BATF","TNFRSF4", "FOXP3","CTLA4","LAIR2"))
		seurat_feature_plot(colon, sample_name, reduction, "T_memory", c("BACH2", "IFNG", "STIM2", "ID2", "IFNAR1", "IL12RB2", "PTPRC"))
		seurat_feature_plot(colon, sample_name, reduction, "T_naive", c("CCR7", "CD28", "ETS1"))
		seurat_feature_plot(colon, sample_name, reduction, "activated_CD4", c("IL4R", "STAT1", "MAL", "SOCS1", "IL2", "ODC1", "WARS"))
		seurat_feature_plot(colon, sample_name, reduction, "T_activated", c("TNF", "IFNG", "JUN", "FOS", "CD69", "REL"))
		seurat_feature_plot(colon, sample_name, reduction, "Th17_CD4", c("IL17A", "CTSH", "KLRB1", "IL26"))
		seurat_feature_plot(colon, sample_name, reduction, "T_exhauseted", c("PDCD1", "HAVCR2", "LAG3","CD101", "CD38", "CXCR6", "TIGIT"))
		seurat_feature_plot(colon, sample_name, reduction, "T_term_exhausted", c("TOX", "GZMB", "ENTPD1", "ITGAE"))
		seurat_feature_plot(colon, sample_name, reduction, "Cycling", c("MKI67", "TOP2A", "CDK1"))
		seurat_feature_plot(colon, sample_name, reduction, "DC", c("ITGAX", "ITGAM", "ITGAE", "LY75"))
	}
	# # Uncomment if you want to save an object and run through azimuth and load the results. 
	# temp <- DietSeurat(colon, assays = c("decontXcounts"))
	# temp <- RenameAssays(object = temp, decontXcounts = 'RNA')
	# saveRDS(temp, "immune_for_azimuth.rds")
	# predictions <- read.delim('./azimuth_pred.tsv', row.names = 1)
	# colon <- AddMetaData(
	# 	object = colon,
	# 	metadata = predictions)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 3) Cluster identification
if (3 %in% execute_steps){
	new.cluster.ids <- c(
		"CD8",
		"CD8",
		"CD4",
		"Mono_Macrophages",
		"Plasma",
		"B Cells",
		"Plasma",
		"NK",
		"B Cells",
		"Plasma",
		"Mono_Macrophages",
		"CD4",
		"Mono_Macrophages",
		"Mono_Macrophages",
		"Mast",
		"T Cells",
		"CyclingImmune",
		"B Cells",
		"CD8",
		"T Cells",
		"Plasma",
		"Mono_Macrophages",
		"ILC",
		"DC",
		"Mono_Macrophages",
		"CD8")

	pal <- c("#faa818")
    names(pal) <- "B Cells"
    pal["CD4"] <- "#fbdf72"
    pal["CD8"] <- "#367d7d"
    pal["CyclingImmune"] <- "#41a30d"
    pal["DC"] <- "#d33502"
    pal["ILC"] <- "#6ebcbc"
    pal["Mono_Macrophages"] <- "#f5b390"
    pal["Mast"] <- "#342739"
    pal["NK"] <- "#bed678"
    pal["Plasma"] <- "#a6d9ee"
    pal["T Cells"] <- "#0d74b6"

	identities <- as.character(colon@meta.data$seurat_clusters)
	for (i in 0:length(new.cluster.ids)){
		identities[identities==as.character(i)] <- new.cluster.ids[i+1]
	}
	colon <- AddMetaData(colon, identities, col.name = "CellType")
	pdf(paste0("./UMAP_cell_type_initial.pdf"), width = 4, height = 5, onefile = F)
	DimPlot(colon, reduction = "umapharmony", pt.size = 0.4, group.by = "CellType", cols = pal) + theme_ArchR()
	dev.off()

	# save
	saveRDS(DietSeurat(
		colon,
		counts = TRUE,
		data = TRUE,
		scale.data = FALSE,
		features = NULL,
		assays = c("RNA", "soupXcounts", "soupXcounts0p2", "soupXcounts0p3", "soupXcounts0p4", "decontXcounts"),
		dimreducs = c("pca", "umap", "harmony", "umapharmony"),
		graphs = NULL
		),"diet_clustered_full_colon_immune_proj_seurat_filtering_complete.rds")
	
	saveRDS(colon, "clustered_full_colon_immune_proj_seurat_filtering_complete.rds")

	colon <- readRDS("clustered_full_colon_immune_proj_seurat_filtering_complete.rds")
}


