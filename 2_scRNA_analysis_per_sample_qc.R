# Script to perform per individual sample RNA qc and filtering for the single cell hubmap colon project.
# WRB 2022

# Import libraries
.libPaths("/oak/stanford/groups/wjg/wbecker/libraries/R_LIBS_4p1p2/")
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

# load helper functions
source("/oak/stanford/groups/wjg/wbecker/other/hubmap_single_cell/scRNA_reproduce/scripts/single_cell_rna_helper_functions.R")

# Define variables
sample_name <- "all_samples"
individual_qc_and_dublet_plot_location <- "/oak/stanford/groups/wjg/wbecker/other/hubmap_single_cell/scRNA_reproduce/all_cells/individual_qc/"
analysis_parent_folder <- "/oak/stanford/groups/wjg/wbecker/other/hubmap_single_cell/scRNA_reproduce/all_cells/"
setwd(analysis_parent_folder)


# set directory for individual processing
setwd("/oak/stanford/groups/wjg/wbecker/other/hubmap_single_cell/scRNA_reproduce/initial_processing/")


# Define sets and locations of files for initial processing
scRNA_data_path <- c(
	"/oak/stanford/groups/wjg/wbecker/scHuBMAP_HTAN/scRNA/cell_ranger_outputs/",
	"/oak/stanford/groups/wjg/wbecker/scHuBMAP_HTAN/scRNA/cell_ranger_outputs/",
	"/oak/stanford/groups/wjg/wbecker/scHuBMAP_HTAN/scRNA/cell_ranger_outputs/",
	"/oak/stanford/groups/wjg/wbecker/scHuBMAP_HTAN/hubmap_multiome_plus_redos/cell_ranger_outputs/rna_redos/",
	"/oak/stanford/groups/wjg/wbecker/scHuBMAP_HTAN/hubmap_multiome_plus_redos/cell_ranger_outputs/")
scRNA_set_names <- c("B001", "B004", "B005", "redos", "B004_multiome")
scRNA_sets <- list(c("B001-A-001","B001-A-006","B001-A-101","B001-A-201", "B001-A-301","B001-A-401","B001-A-406", "B001-A-501"),
			c("B004-A-004", "B004-A-008", "B004-A-104", "B004-A-204", "B004-A-304","B004-A-504","B004-A-404","B004-A-408"),#,"B004-A-404","B004-A-408"
			c("B005-A-001","B005-A-002","B005-A-101","B005-A-201","B005-A-301","B005-A-401","B005-A-402","B005-A-501"), #,"B005-A-501"
			c("B004-A-104-R2", "B004-A-408-R2", "B004-A-504-R2", "B005-A-501-R2"),
			c("B004-A-404-R2"))

multiome_data_path <- c(
	"/oak/stanford/groups/wjg/wbecker/scHuBMAP_HTAN/hubmap_multiome_plus_redos/cell_ranger_outputs/",
	"/oak/stanford/groups/wjg/wbecker/scHuBMAP_HTAN/hubmap_multiome_plus_redos/cell_ranger_outputs/",
	"/oak/stanford/groups/wjg/wbecker/scHuBMAP_HTAN/hubmap_multiome_plus_redos/cell_ranger_outputs/",
	"/oak/stanford/groups/wjg/wbecker/scHuBMAP_HTAN/hubmap_multiome_plus_redos/cell_ranger_outputs/",
	"/oak/stanford/groups/wjg/wbecker/scHuBMAP_HTAN/hubmap_multiome_plus_redos/cell_ranger_outputs/",
	"/oak/stanford/groups/wjg/wbecker/scHuBMAP_HTAN/hubmap_multiome_plus_redos/cell_ranger_outputs/")
multiome_set_names <- c("B006", "S1", "S2", "S3", "S4", "S5")
multiome_sets <- list(c("B006-A-002", "B006-A-201", "B006-A-401", "B006-A-402","B006-A-201-R2", "B006-A-101", "B006-A-001", "B006-A-501", "B006-A-301"),
	c("Bei_CO1-1_B010-A-101", "Bei_CO1-2_B010-A-201", "Bei_CO1-3_B010-A-405", "Bei_CO1-4_B010-A-501", "Bei_CO1-5_B010-A-001","Bei_CO1-6_B010-A-002", "Bei_CO1-7_B010-A-301", "Bei_CO1-8_B010-A-401"), 
	c("Bei_CO2-1_B011-A-001", "Bei_CO2-10_B012-A-002","Bei_CO2-11_B012-A-301", "Bei_CO2-12_B012-A-401", "Bei_CO2-13_B012-A-101", "Bei_CO2-14_B012-A-201","Bei_CO2-17_B008-A-002", "Bei_CO2-18_B008-A-101"),
	c("Bei_CO2-19_B008-A-201", "Bei_CO2-2_B011-A-002", "Bei_CO2-20_B008-A-401", "Bei_CO2-21_B008-A-402", "Bei_CO2-22_B008-A-501", "Bei_CO2-23_B012-A-405", "Bei_CO2-24_B012-A-501", "Bei_CO2-25_B009-A-001"), 
	c("Bei_CO2-26_B008-A-001", "Bei_CO2-27_B009-A-301", "Bei_CO2-28_B009-A-401", "Bei_CO2-29_B009-A-101", "Bei_CO2-3_B011-A-301", "Bei_CO2-30_B008-A-301", "Bei_CO2-31_B009-A-405"), #, "Bei_CO2-32_B009-A-501"
	c("Bei_CO2-4_B011-A-401", "Bei_CO2-5_B011-A-101","Bei_CO2-6_B011-A-201", "Bei_CO2-7_B011-A-405", "Bei_CO2-8_B011-A-501", "Bei_CO2-9_B012-A-001"))

# load atac cells passing filters
atac_cells <- read.table("/oak/stanford/groups/wjg/wbecker/other/hubmap_single_cell/scATAC/projects/initial_post_filter_atac_cells.txt")

# load immune cell types for these samples (from downstream analysis) if desired for examining how well ambient RNA removal is working
celltypes <- read.table("/oak/stanford/groups/wjg/wbecker/other/hubmap_single_cell/scRNA/v1/immune/counts_norm_scale_harmony/immune_cell_types.tsv")

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# 1) Create seurat objects for individual 10x runs, run doblet finder and filter most likely doublets, merge into a seurat object containing all samples
# First process the singleome rna-seq samples

# The stepwise merges into smaller groups were done because we ran into memory issues with merging all of the samples at once.

for (j in 1:length(scRNA_set_names)){
	samples <- scRNA_sets[[j]]
	print(paste0(scRNA_data_path[j], samples[1], "/outs/raw_feature_bc_matrix/"))
	data_directory <- paste0(scRNA_data_path[j], samples[1], "/outs/raw_feature_bc_matrix/")
	sample1 <- make_seurat_object_and_doublet_removal(data_directory, samples[1], runDoubletFinder = TRUE, runSoupX = TRUE, cell_filter = NULL, passed_immune_labels = celltypes)
	seu_list <- c()
	if (length(samples)>1){
		for (i in 2:length(samples)){
			data_directory <- paste0(scRNA_data_path[j], samples[i], "/outs/raw_feature_bc_matrix/")
			seu_list <- c(seu_list, make_seurat_object_and_doublet_removal(data_directory, samples[i], runDoubletFinder = TRUE, runSoupX = TRUE, cell_filter = NULL, passed_immune_labels = celltypes))
			while (!is.null(dev.list()))  dev.off()
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
	print(paste0(multiome_data_path[j], samples[1], "/outs/raw_feature_bc_matrix/"))
	data_directory <- paste0(multiome_data_path[j], samples[1], "/outs/raw_feature_bc_matrix/")
	# get cells from atac qc and doublet filtering
	sample_name_only <- rev(str_split(samples[1], "_")[[1]])[1]
	current_cells <- str_replace(atac_cells$x[grepl(sample_name_only, atac_cells$x)], paste0(sample_name_only, "#"), "")
	sample1 <- make_seurat_object_and_doublet_removal(data_directory, sample_name_only, runDoubletFinder = FALSE, runSoupX = TRUE, cell_filter = current_cells)
	seu_list <- c()
	if (length(samples)>1){
		for (i in 2:length(samples)){
			data_directory <- paste0(multiome_data_path[j], samples[i], "/outs/raw_feature_bc_matrix/")
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
# Fix cell names that got messed up and combine into a single project
colonP1 <- readRDS("non_multiome_colon_proj_seurat.rds")
colon <- readRDS("multiome_colon_proj_seurat.rds")

# add missing rowname info
P1_names <- Cells(colonP1)
P1_names[substr(P1_names, 1, nchar(P1_names)-19) == ""] <- paste0("B004-A-404-R2_", P1_names[substr(P1_names, 1, nchar(P1_names)-19) == ""])
colonP1 <- RenameCells(colonP1,new.names = P1_names)

P2_names <- Cells(colon)
P2_names <- substr(P2_names, nchar(P2_names)-28, nchar(P2_names))
P2_names[substr(P2_names, 1, nchar(P2_names)-19) == "6-A-201-R2"] <- paste0("B00", P2_names[substr(P2_names, 1, nchar(P2_names)-19) == "6-A-201-R2"])

colonP1 <- RenameCells(colonP1,new.names = P1_names)
colon <- RenameCells(colon,new.names = P2_names)

colon <- merge(colonP1, y = colon, project = "full_colon_project")
saveRDS(colon, "full_intestine_proj_seurat.rds")

