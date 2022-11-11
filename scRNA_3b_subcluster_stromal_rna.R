# Script to analyze stromal scRNA cells from HuBMAP intestine data
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
sample_name <- "all_samples" # used as label for saving plots
analysis_parent_folder <- paste0(parent_folder, "scRNA_reproduce/stromal/decontx_norm_scale_harmony/") # folder for this analysis
setwd(analysis_parent_folder)

# load helper functions
source(paste0(parent_folder, "scRNA_reproduce/scripts/single_cell_rna_helper_functions.R"))

# Set steps to run
execute_steps <- c(1,2,3,4,5,6,7)

# Load previously defined seurat object that contains only the cells classified as immune cells
colon <- readRDS(paste0(parent_folder, "scRNA_reproduce/all_cells/decontx_norm_scale_harmony/colon_stromal_all_samples_initial.rds"))
extra_cells <- readRDS(paste0(parent_folder, "scRNA_reproduce/immune/decontx_norm_scale_harmony/immune_inital_move_to_stromal.rds"))
colon <- merge(colon, y = c(extra_cells), project = "full_stromal")

# Markers from Regev Lab paper
RNA_markers_epithelial <- read.csv(file = paste0(parent_folder, "RegevLabMarkers/S2_Regev_Cell_epithelial.csv"))
RNA_markers_immune <- read.csv(file = paste0(parent_folder, "RegevLabMarkers/S2_Regev_Cell_immune.csv"))
RNA_markers_stromal <- read.csv(file = paste0(parent_folder, "RegevLabMarkers/S2_Regev_Cell_stromal.csv"))
RNA_markers <- rbind(RNA_markers_epithelial, RNA_markers_immune, RNA_markers_stromal)

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
if (1 %in% execute_steps){
	colon <- seurat_standard_normalize_and_scale(colon, cluster = FALSE, ndims <- 20)
	plotUMAP(colon, reduction = "umap", save_name = "pre_filter")
	colon <- run_harmony(colon)
	set.seed(1)
	colon <- FindClusters(colon, resolution = 0.7)
	plotUMAP(colon, reduction = "umapharmony", save_name = "pre_filter")
	
	# ID Cluster Markers and look at overlapping specific markers and print out overlapping cases
	colon.markers <- FindAllMarkers(colon, only.pos = TRUE, min.pct = 0.05, max.cells.per.ident = 250)
	RNA_markers_spec <- RNA_markers[RNA_markers$either,]
	nclusters <- length(unique(colon.markers$cluster))-1
	for (val in 0:nclusters){
	    clusterSpecficMarkers <- RNA_markers_spec[RNA_markers_spec$gene %in% colon.markers[colon.markers$cluster == val,]$gene,]
	    print(val)
	    print(clusterSpecficMarkers[,c("ident", "gene")])
	}
	RNA_markers_spec[RNA_markers_spec$gene %in% colon.markers[colon.markers$cluster == val,]$gene,c("ident")]
	
	# Label Transfer and singleR
	colon <- run_singleR_and_plot(colon, types_to_use = c("Smooth_muscle_cells","Endothelial_cells","Neurons","Fibroblasts","Neuroepithelial_cell","Epithelial_cells","T_cells","B_cells"))

	# Plot marker genes
	sample_name_temp <- paste0(sample_name, "_pre_fiter")
	reductions_to_plot <- c("umapharmony")
	for (reduction in reductions_to_plot){
		seurat_feature_plot(colon, sample_name_temp, reduction, "Fibroblasts", c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "FAP", "CBLN2", "SPOCK1", "ACSS3"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Fibroblasts_cell_subtypes", c("RSPO3", "CCL11", "WNT5B", "BMP4", "CHI3L1", "ACTA2", "WNT2B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Myofibroblasts", c("SYT10", "SOSTDC1", "DES", "MYH11", "TAGLN", "ACTA2", "TPM4"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "stromal_other", c("FAM110D", "INHBB", "NPR1", "NOVA2", "GPIHBP1", "SOX17", "VWF", "PLVAP", "CDH5", "S100B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Inflammatory_fibroblasts", c("CHI3L1","MMP3","PLAU","MMP1","TRAFD1","GBP1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Post_capillary_venules", c("SELP", "ZNF385D", "FAM155A", "GALNT15", "MADCAM1", "CORT"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Lymph_endothelial", c("PROX1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Pericytes", c("MCAM", "COX4I2", "KCNJ8", "HIGD1B", "RGS5", "NOTCH3", "HEYL", "FAM162B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Microvascular", c("PLVAP","CD36","DYSF","NRP1","SH3BP5","EXOC3L2","FABP5","VWA1","BAALC","PRSS23","RAPGEF4","APLN","HTRA1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "SchwanCell", c("S100A1", "SOX10", "EGR2", "MBP", "MPZ", "GAP43", "NCAM", "P75NTR"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Nerve", c("MAP2", "RBFOX3", "DLG4", "SYP"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Stem_Maintenance", c("CD34", "GREM1", "RSPO3", "WNT2B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Epithelial_Growth", c("SEMA3B", "SEMA3E", "SEMA3C"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "BMP_Signaling", c("BMP5", "BMP4", "BMP2"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Immune_attraction", c("CXCL12", "CXCL14"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Eos_recruitment", c("CCL11", "CCL13"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "myeloid_recruitment", c("CCL2", "CCL8"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "inflammatory_CAF_markers", c("TWIST1", "WNT2", "FAP"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "ECM", c("FBLN1", "PCOLCE2", "MFAP5"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "basement_membrane", c("COL4A5", "COL4A6"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "PITX", c("PITX1", "PITX2", "PITX3"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Insterstitial_Cells_Cajal", c("KIT", "ANO1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Myofib_sub", c("TAGLN", "ACTA2", "DCN", "DES", "HHIP", "NPNT", "SYT10", "RSPO2", "SYT1", "PTGER1", "CNN1"))
	}
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 2) Remove possible doublet clusters and redo analysis x2.
if (2 %in% execute_steps){
	# We will move possible doublet clusters
	# Save ids of immune, epithelial, and stromal cells
	immune_clusters <- c(8,13) 
	epithelial_clusters <- c(1)

	# Subset to make immune, stromal, and epithelial projects
	colon_epitehlial <- DietSeurat(subset(colon, subset = seurat_clusters %in% epithelial_clusters))
	saveRDS(colon_epitehlial, file = "stromal_inital_move_to_epithelial.rds")

	colon_immune <- DietSeurat(subset(colon, subset = seurat_clusters %in% immune_clusters))
	saveRDS(colon_immune, file = "stromal_inital_likely_stromal_immune_doublets.rds")

	sum(colnames(colon_immune) %in% colnames(extra_cells))

	bad_clusters <- c(1,8,13)

	# Subset
	colon_new <- DietSeurat(subset(colon, subset = seurat_clusters %ni% bad_clusters))

	# Redo the initial analysis steps
	colon_new <- seurat_standard_normalize_and_scale(colon_new, cluster = FALSE, ndims <- 20)
	plotUMAP(colon_new, reduction = "umap", save_name = "post_filter_1")
	colon <- colon_new
	set.seed(1)
	colon <- run_harmony(colon, resolution = 1.0)
	plotUMAP(colon, reduction = "umapharmony", save_name = "post_filter_1")
	saveRDS(colon, file = "stromal_all_samples_clustered_filtered_post_filter1.rds")

	sample_name_temp <- paste0(sample_name, "_post_doublet_removal")
	reductions_to_plot <- c("umapharmony")
	for (reduction in reductions_to_plot){
		seurat_feature_plot(colon, sample_name_temp, reduction, "Fibroblasts", c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "FAP", "CBLN2", "SPOCK1", "ACSS3"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Fibroblasts_cell_subtypes", c("RSPO3", "CCL11", "WNT5B", "BMP4", "CHI3L1", "ACTA2", "WNT2B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Myofibroblasts", c("SYT10", "SOSTDC1", "DES", "MYH11", "TAGLN", "ACTA2", "TPM4"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "stromal_other", c("FAM110D", "INHBB", "NPR1", "NOVA2", "GPIHBP1", "SOX17", "VWF", "PLVAP", "CDH5", "S100B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Inflammatory_fibroblasts", c("CHI3L1","MMP3","PLAU","MMP1","TRAFD1","GBP1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Post_capillary_venules", c("SELP", "ZNF385D", "FAM155A", "GALNT15", "MADCAM1", "CORT"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Lymph_endothelial", c("PROX1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Pericytes", c("MCAM", "COX4I2", "KCNJ8", "HIGD1B", "RGS5", "NOTCH3", "HEYL", "FAM162B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Microvascular", c("PLVAP","CD36","DYSF","NRP1","SH3BP5","EXOC3L2","FABP5","VWA1","BAALC","PRSS23","RAPGEF4","APLN","HTRA1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "SchwanCell", c("S100A1", "SOX10", "EGR2", "MBP", "MPZ", "GAP43", "NCAM", "P75NTR"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Nerve", c("MAP2", "RBFOX3", "DLG4", "SYP"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Stem_Maintenance", c("CD34", "GREM1", "RSPO3", "WNT2B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Epithelial_Growth", c("SEMA3B", "SEMA3E", "SEMA3C"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "BMP_Signaling", c("BMP5", "BMP4", "BMP2"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Immune_attraction", c("CXCL12", "CXCL14"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Eos_recruitment", c("CCL11", "CCL13"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "myeloid_recruitment", c("CCL2", "CCL8"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "inflammatory_CAF_markers", c("TWIST1", "WNT2", "FAP"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "ECM", c("FBLN1", "PCOLCE2", "MFAP5"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "basement_membrane", c("COL4A5", "COL4A6"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "PITX", c("PITX1", "PITX2", "PITX3"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Insterstitial_Cells_Cajal", c("KIT", "ANO1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Myofib_sub", c("TAGLN", "ACTA2", "DCN", "DES", "HHIP", "NPNT", "SYT10", "RSPO2", "SYT1", "PTGER1", "CNN1"))
	}

	colon.markers <- FindAllMarkers(colon, only.pos = TRUE, min.pct = 0.05, max.cells.per.ident = 250)
}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 3) Clean up a few more low quality clusters
if (3 %in% execute_steps){
	# Save ids of immune, epithelial, and stromal cells
	immune_clusters <- c() 
	epithelial_clusters <- c(25)

	# Subset to make immune, stromal, and epithelial projects
	colon_epitehlial <- DietSeurat(subset(colon, subset = seurat_clusters %in% epithelial_clusters))
	saveRDS(colon_epitehlial, file = "stromal_inital_move_to_epithelial_R2.rds")

	bad_clusters <- c(immune_clusters, epithelial_clusters)
	colon_new <- DietSeurat(subset(colon, subset = seurat_clusters %ni% bad_clusters))

	# Subset
	colon_new <- seurat_standard_normalize_and_scale(colon_new, cluster = FALSE, ndims <- 20)
	plotUMAP(colon_new, reduction = "umap", save_name = "post_filter_2")
	colon <- colon_new
	set.seed(1)
	colon <- run_harmony(colon, resolution = 1.0)
	set.seed(1)
	colon <- FindClusters(colon, resolution = 1.0)
	plotUMAP(colon, reduction = "umapharmony", save_name = "post_filter_2")

	sample_name_temp <- paste0(sample_name, "_post_doublet_removal_2")
	reductions_to_plot <- c("umapharmony")
	for (reduction in reductions_to_plot){
		seurat_feature_plot(colon, sample_name_temp, reduction, "Fibroblasts", c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "FAP", "CBLN2", "SPOCK1", "ACSS3"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Fibroblasts_cell_subtypes", c("RSPO3", "CCL11", "WNT5B", "BMP4", "CHI3L1", "ACTA2", "WNT2B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Myofibroblasts", c("SYT10", "SOSTDC1", "DES", "MYH11", "TAGLN", "ACTA2", "TPM4"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "stromal_other", c("FAM110D", "INHBB", "NPR1", "NOVA2", "GPIHBP1", "SOX17", "VWF", "PLVAP", "CDH5", "S100B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Inflammatory_fibroblasts", c("CHI3L1","MMP3","PLAU","MMP1","TRAFD1","GBP1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Post_capillary_venules", c("SELP", "ZNF385D", "FAM155A", "GALNT15", "MADCAM1", "CORT"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Lymph_endothelial", c("PROX1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Pericytes", c("MCAM", "COX4I2", "KCNJ8", "HIGD1B", "RGS5", "NOTCH3", "HEYL", "FAM162B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Microvascular", c("PLVAP","CD36","DYSF","NRP1","SH3BP5","EXOC3L2","FABP5","VWA1","BAALC","PRSS23","RAPGEF4","APLN","HTRA1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "SchwanCell", c("S100A1", "SOX10", "EGR2", "MBP", "MPZ", "GAP43", "NCAM", "P75NTR"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Nerve", c("MAP2", "RBFOX3", "DLG4", "SYP", "RBFOX1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Stem_Maintenance", c("CD34", "GREM1", "RSPO3", "WNT2B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Epithelial_Growth", c("SEMA3B", "SEMA3E", "SEMA3C"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "BMP_Signaling", c("BMP5", "BMP4", "BMP2"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Immune_attraction", c("CXCL12", "CXCL14"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Eos_recruitment", c("CCL11", "CCL13"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "myeloid_recruitment", c("CCL2", "CCL8"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "inflammatory_CAF_markers", c("TWIST1", "WNT2", "FAP"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "ECM", c("FBLN1", "PCOLCE2", "MFAP5"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "basement_membrane", c("COL4A5", "COL4A6"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Insterstitial_Cells_Cajal", c("KIT", "ANO1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "PITX", c("PITX1", "PITX2", "PITX3"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Myofib_sub", c("TAGLN", "ACTA2", "DCN", "DES", "HHIP", "NPNT", "SYT10", "RSPO2", "SYT1", "PTGER1", "CNN1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "cycling", c("MKI67", "TOP2A", "CDK1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Adipocytes", c("LPL", "PLIN1"))
	}


	# Save ids of immune, epithelial, and stromal cells
	immune_clusters <- c() 
	epithelial_clusters <- c(25)

	# Subset to make immune, stromal, and epithelial projects
	colon_epitehlial <- DietSeurat(subset(colon, subset = seurat_clusters %in% epithelial_clusters))
	saveRDS(colon_epitehlial, file = "stromal_inital_move_to_epithelial_R3.rds")

	bad_clusters <- c(immune_clusters, epithelial_clusters)
	colon_new <- DietSeurat(subset(colon, subset = seurat_clusters %ni% bad_clusters))

	# Subset
	colon_new <- seurat_standard_normalize_and_scale(colon_new, cluster = FALSE, ndims <- 20)
	plotUMAP(colon_new, reduction = "umap", save_name = "post_filter_3")
	colon <- colon_new
	set.seed(1)
	colon <- run_harmony(colon, resolution = 1.0)
	set.seed(1)
	colon <- FindClusters(colon, resolution = 0.9)
	plotUMAP(colon, reduction = "umapharmony", save_name = "post_filter_3")

	sample_name_temp <- paste0(sample_name, "_post_doublet_removal_3")
	reductions_to_plot <- c("umapharmony")
	for (reduction in reductions_to_plot){
		seurat_feature_plot(colon, sample_name_temp, reduction, "Fibroblasts", c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "FAP", "CBLN2", "SPOCK1", "ACSS3"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Fibroblasts_cell_subtypes", c("RSPO3", "CCL11", "WNT5B", "BMP4", "CHI3L1", "ACTA2", "WNT2B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Myofibroblasts", c("SYT10", "SOSTDC1", "DES", "MYH11", "TAGLN", "ACTA2", "TPM4"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "stromal_other", c("FAM110D", "INHBB", "NPR1", "NOVA2", "GPIHBP1", "SOX17", "VWF", "PLVAP", "CDH5", "S100B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Inflammatory_fibroblasts", c("CHI3L1","MMP3","PLAU","MMP1","TRAFD1","GBP1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Post_capillary_venules", c("SELP", "ZNF385D", "FAM155A", "GALNT15", "MADCAM1", "CORT"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Lymph_endothelial", c("PROX1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Pericytes", c("MCAM", "COX4I2", "KCNJ8", "HIGD1B", "RGS5", "NOTCH3", "HEYL", "FAM162B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Microvascular", c("PLVAP","CD36","DYSF","NRP1","SH3BP5","EXOC3L2","FABP5","VWA1","BAALC","PRSS23","RAPGEF4","APLN","HTRA1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "SchwanCell", c("S100A1", "SOX10", "EGR2", "MBP", "MPZ", "GAP43", "NCAM", "P75NTR"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Nerve", c("MAP2", "RBFOX3", "DLG4", "SYP", "RBFOX1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Stem_Maintenance", c("CD34", "GREM1", "RSPO3", "WNT2B"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Epithelial_Growth", c("SEMA3B", "SEMA3E", "SEMA3C"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "BMP_Signaling", c("BMP5", "BMP4", "BMP2"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Immune_attraction", c("CXCL12", "CXCL14"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Eos_recruitment", c("CCL11", "CCL13"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "myeloid_recruitment", c("CCL2", "CCL8"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "inflammatory_CAF_markers", c("TWIST1", "WNT2", "FAP"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "ECM", c("FBLN1", "PCOLCE2", "MFAP5"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "basement_membrane", c("COL4A5", "COL4A6"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Insterstitial_Cells_Cajal", c("KIT", "ANO1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "PITX", c("PITX1", "PITX2", "PITX3"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Myofib_sub", c("TAGLN", "ACTA2", "DCN", "DES", "HHIP", "NPNT", "SYT10", "RSPO2", "SYT1", "PTGER1", "CNN1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "cycling", c("MKI67", "TOP2A", "CDK1"))
		seurat_feature_plot(colon, sample_name_temp, reduction, "Adipocytes", c("LPL", "PLIN1"))
	}

	pdf(paste0("./post_doub_rem_UMAP_predicted_id_label_transfer-hubmapstromal.pdf"), width = 12)
	print(DimPlot(colon, reduction = "umapharmony", group.by = "predicted.id.hubmapstromal", cols = paletteDiscrete(values = unique(colon@meta.data$predicted.id.hubmapstromal), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./post_doub_rem_UMAP_singleR.pdf"), width = 12)
	print(DimPlot(colon, reduction = "umapharmony", group.by = "SingleR.labels", cols = paletteDiscrete(values = unique(colon@meta.data$SingleR.labels), set = "stallion", reverse = FALSE)))
	dev.off()

	colon <- run_singleR_and_plot(colon, types_to_use = c("Smooth_muscle_cells","Endothelial_cells","Neurons","Fibroblasts","Neuroepithelial_cell","Epithelial_cells","T_cells","B_cells"))

}

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 4) Initial cluster identification
if (4 %in% execute_steps){
	new.cluster.ids <- c(
		"Myofibroblasts/SM 1",
		"Myofibroblasts/SM 1",
		"Crypt Fibroblasts 1 WNT2B+",
		"Crypt Fibroblasts 3 RSPO3+",
		"Myofibroblasts/SM DES High",
		"Crypt Fibroblasts 2",
		"Endothelial-CD36+ Microvascular",
		"Myofibroblasts/SM DES High",
		"Myofibroblasts/SM 3",
		"Glia",
		"Villus Fibroblasts WNT5B+",
		"Endothelial-Venules",
		"Lymphatic Endothelial Cells",
		"Pericytes",
		"Myofibroblasts/SM 2",
		"Myofibroblasts/SM 1",
		"Myofibroblasts/SM 1",
		"Crypt Fibroblasts 1 WNT2B+",
		"Crypt Fibroblasts 1 WNT2B+",
		"ICC",
		"Crypt Fibroblasts 3 RSPO3+",
		"Endothelial-CD36+ Microvascular",
		"Myofibroblasts/SM 1",
		"Endothelial-CD36+ Microvascular",
		"Neurons",
		"Crypt Fibroblasts 3 RSPO3+",
		"Adipocytes",
		"Crypt Fibroblasts 1 WNT2B+",
		"Crypt Fibroblasts 3 RSPO3+")
	
	identities <- as.character(colon@meta.data$seurat_clusters)
	for (i in 0:length(new.cluster.ids)){
		identities[identities==as.character(i)] <- new.cluster.ids[i+1]
	}
	colon <- AddMetaData(colon, identities, col.name = "CellType")

	pal <- c()
    pal["Adipocytes"] <- "#93AA00"
    pal["Crypt Fibroblasts 1 WNT2B+"] <- "#FFB300"
    pal["Crypt Fibroblasts 2"] <- "#803E75"
    pal["Crypt Fibroblasts 3 RSPO3+"] <- "#FF6800"
    #pal["Crypt Fibroblasts 4 RSPO3+"] <- "#F13A13"
    pal["Endothelial-CD36+ Microvascular"] <- "#C10020"
    pal["Endothelial-Venules"] <- "#817066"
    pal["Lymphatic Endothelial Cells"] <- "#A6BDD7"
    pal["Glia"] <- "#CEA262"
    pal["Myofibroblasts/SM 1"] <- "#007D34"
    pal["Myofibroblasts/SM 2"] <- "#F6768E"
    pal["Myofibroblasts/SM 3"] <- "#00538A"
    #pal["Myofibroblasts/SM 4"] <- "#FF7A5C"
    pal["Myofibroblasts/SM DES High"] <- "#FF8E00"
    pal["Neurons"] <- "#F4C800"
    pal["Pericytes"] <- "#53377A"
    pal["ICC"] <- "#7F180D"
    pal["Villus Fibroblasts WNT5B+"] <- "#B32851"

	pdf(paste0("./UMAP_cell_type_vector.pdf"), width = 6.5, onefile = F)
	DimPlot(colon, reduction = "umapharmony", group.by = "CellType", cols = pal)+theme_ArchR()
	dev.off()

    # save
	saveRDS(colon, file = "stromal_all_samples_clustered_filtered.rds")
	saveRDS(DietSeurat(
		colon,
		counts = TRUE,
		data = TRUE,
		scale.data = FALSE,
		features = NULL,
		assays = NULL,
		dimreducs = c("pca", "umap", "harmony", "umapharmony"),
		graphs = NULL
		),"diet_stromal_all_samples_clustered_filtered.rds")


	# Make some fibroblast specific dotplots
	fib_only <- colon
	Idents(fib_only) <- "CellType"

	fib_only@active.ident <- factor(fib_only@active.ident, 
	                            levels=c("Adipocytes",
	                            	"Neurons", "Glia", 
	                            	"ICC",
	                            	"Lymphatic Endothelial Cells", "Endothelial-CD36+ Microvascular", "Endothelial-Venules", 
	                            	"Pericytes", "Myofibroblasts/SM DES High", "Myofibroblasts/SM 1","Myofibroblasts/SM 2", "Myofibroblasts/SM 3", #"Myofibroblasts/SM 4",
	                            	"Villus Fibroblasts WNT5B+", 
	                            	"Crypt Fibroblasts 1 WNT2B+", "Crypt Fibroblasts 2", "Crypt Fibroblasts 3 RSPO3+"#, "Crypt Fibroblasts 4 RSPO3+"
	                            ))

	markers <- c(
		"PLIN1", "LPL", # Adipocytes
		"SYP", "SYT1", "MIR124-2HG", "RBFOX1", # Neuron
		"SOX10", "CDH19", "PLP1", "S100B", # Enteric Glia
		"KIT", "ANO1", #ICC
		"PROX1", "RELN", # Lymph Endothelial
		"CD36", "PECAM1", #Endothelial
		"SELP", "MADCAM1",
		"NOTCH3", "MCAM", "RGS5", #Pericytes
		"TAGLN", "ACTA2", "MYH11", # Myofibroblasts , "TPM4"
		"DES", #Smooth Muscle
		"COL6A1", "COL6A2", # General Fib Markers
		"WNT5B", # Villus WNT Signaling (WNT antagonist)
		"BMP5", "BMP4", "BMP2", "PDGFRA", # BMP signaling
		"CCL8", "CCL13", "CCL19", # Immune recruitment "CCL21", 
		"ADAMDEC1", 
		"CCL11", "KCNN3", 
		"WNT2B", "RSPO3", "GREM1",  # Stem maintence
		"SEMA3E", #Epithelial growth "SEMA3C"
		"C3"
		#"DCN", "HHIP", "NPNT", "RSPO2", "PTGER1", "CNN1",
		)
	p <- DotPlot(fib_only, features = markers, col.max = 2, dot.min = 0.02) + RotatedAxis()+scale_colour_gradientn(colours = wes_palette("Zissou1", 10, type = "continuous"))
	ggsave(paste0("./", "dot_plot_", sample_name, "_", "Fibroblasts_wo_COL4A_height6" ,".pdf"), plot = p, width = 18, height = 6, useDingbats=FALSE)
}







