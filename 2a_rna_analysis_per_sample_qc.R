# Script to perform per individual sample RNA qc and filtering for the single cell hubmap colon project.
# WRB 2022

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
#library(celda)
`%notin%` <- Negate(`%in%`)
set.seed(1)

# Define variables
parent_directory <- "/hubmap_single_cell/"
sample_name <- "all_samples"
individual_qc_and_dublet_plot_location <- paste0(parent_directory, "scRNA/all_cells/individual_qc/")
analysis_parent_folder <- paste0(parent_directory, "scRNA/initial_processing/")
setwd(analysis_parent_folder)

# load helper functions
source(paste0(parent_directory, "scRNA/scripts/single_cell_rna_helper_functions.R"))

# load atac cells passing filters
atac_cells <- read.table("scATAC/projects/initial_post_filter_atac_cells.txt")

# load immune cell types for these samples (from downstream analysis) if desired for examining how well ambient RNA removal is working
celltypes <- NULL

# Define sets and locations of files for initial processing
scRNA_data_path <- paste0(parent_directory, "data/scRNA/cell_ranger_outputs/")
scRNA_set_names <- c("B001", "B004", "B005", "redos")
scRNA_sets <- list(c("B001-A-001","B001-A-006","B001-A-101","B001-A-201", "B001-A-301","B001-A-401","B001-A-406", "B001-A-501"),
			c("B004-A-004", "B004-A-008", "B004-A-104", "B004-A-204", "B004-A-304","B004-A-504","B004-A-404","B004-A-408"),
			c("B005-A-001","B005-A-002","B005-A-101","B005-A-201","B005-A-301","B005-A-401","B005-A-402","B005-A-501"),
			c("B004-A-104-R2", "B004-A-408-R2", "B004-A-504-R2", "B005-A-501-R2", "B004-A-404-R2"))

multiome_data_path <- paste0(parent_directory, "data/hubmap_multiome_plus_redos/cell_ranger_outputs/")
multiome_set_names <- c("B006", "S1", "S2", "S3", "S4", "S5")
multiome_sets <- list(c("B006-A-002", "B006-A-201", "B006-A-401", "B006-A-402","B006-A-201-R2", "B006-A-101", "B006-A-001", "B006-A-501", "B006-A-301"),
	c("B010-A-101", "B010-A-201", "B010-A-405", "B010-A-501", "B010-A-001", "B010-A-002", "B010-A-301", "B010-A-401"), 
	c("B011-A-001", "B012-A-002", "B012-A-301", "B012-A-401", "B012-A-101", "B012-A-201", "B008-A-002", "B008-A-101"),
	c("B008-A-201", "B011-A-002", "B008-A-401", "B008-A-402", "B008-A-501", "B012-A-405", "B012-A-501", "B009-A-001"), 
	c("B008-A-001", "B009-A-301", "B009-A-401", "B009-A-101", "B011-A-301", "B008-A-301", "B009-A-405"),
	c("B011-A-401", "B011-A-101", "B011-A-201", "B011-A-405", "B011-A-501", "B012-A-001"))

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# Create seurat objects for individual 10x runs, run doblet finder and filter most likely doublets, merge into a seurat object containing all samples
# note that that the merging is done in small batches and then the batches are merged together because I had fewer memory problems with this compared to just merging everything together
# First process the singleome rna-seq samples
for (j in 1:length(scRNA_set_names)){
	samples <- scRNA_sets[[j]]
	print(paste0(scRNA_data_path, samples[1], "/outs/raw_feature_bc_matrix/"))
	data_directory <- paste0(scRNA_data_path, samples[1], "/outs/raw_feature_bc_matrix/")
	sample1 <- make_seurat_object_and_doublet_removal(data_directory, samples[1], runDoubletFinder = TRUE, runSoupX = TRUE, cell_filter = NULL, passed_immune_labels = celltypes)
	seu_list <- c()
	if (length(samples)>1){
		for (i in 2:length(samples)){
			data_directory <- paste0(scRNA_data_path, samples[i], "/outs/raw_feature_bc_matrix/")
			seu_list <- c(seu_list, make_seurat_object_and_doublet_removal(data_directory, samples[i], runDoubletFinder = TRUE, runSoupX = TRUE, cell_filter = NULL, passed_immune_labels = celltypes))
		}
		current_merge <- merge(sample1, y = seu_list, add.cell.ids = samples, project = scRNA_set_names[j])
	} else {
		current_merge <- sample1
	}
	if (j==1){
		colon <- current_merge
	} else if (j>1){
		colon <- merge(colon, y = current_merge, project = "full_colon_project")
	}
}
saveRDS(colon, "non_multiome_colon_proj_seurat.rds")
colon_non_mult <- colon

# Process the multiome samples using the doublet and qc info from the atac data
for (j in 1:length(multiome_set_names)){
	samples <- multiome_sets[[j]]
	print(paste0(multiome_data_path, samples[1], "/outs/raw_feature_bc_matrix/"))
	data_directory <- paste0(multiome_data_path, samples[1], "/outs/raw_feature_bc_matrix/")
	# get cells from atac qc and doublet filtering
	sample_name_only <- rev(str_split(samples[1], "_")[[1]])[1]
	current_cells <- str_replace(atac_cells$x[grepl(sample_name_only, atac_cells$x)], paste0(sample_name_only, "#"), "")
	sample1 <- make_seurat_object_and_doublet_removal(data_directory, sample_name_only, runDoubletFinder = FALSE, runSoupX = TRUE, cell_filter = current_cells)
	seu_list <- c()
	if (length(samples)>1){
		for (i in 2:length(samples)){
			data_directory <- paste0(multiome_data_path, samples[i], "/outs/raw_feature_bc_matrix/")
			# get cells from atac qc and doublet filtering
			sample_name_only <- rev(str_split(samples[i], "_")[[1]])[1]
			current_cells <- str_replace(atac_cells$x[grepl(sample_name_only, atac_cells$x)], paste0(sample_name_only, "#"), "")
			if (length(current_cells) == 0){
				seu_list <- c(seu_list, make_seurat_object_and_doublet_removal(data_directory, sample_name_only, runDoubletFinder = TRUE, runSoupX = TRUE, cell_filter = NULL, passed_immune_labels = celltypes))
				} else {
				seu_list <- c(seu_list, make_seurat_object_and_doublet_removal(data_directory, sample_name_only, runDoubletFinder = FALSE, runSoupX = TRUE, cell_filter = current_cells, passed_immune_labels = celltypes))
			}
		}
		current_merge <- merge(sample1, y = seu_list, add.cell.ids = samples, project = multiome_set_names[j])
	} else {
		current_merge <- sample1
	}
	if (j==1){
		colon <- current_merge
	} else if (j>1){
		colon <- merge(colon, y = current_merge, project = "full_colon_project")
	}
}
saveRDS(colon, "multiome_colon_proj_seurat.rds")

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# Combine into a single project
colon <- merge(colon_non_mult, y = colon, project = "full_colon_project")
saveRDS(colon, "full_intestine_proj_seurat.rds")


