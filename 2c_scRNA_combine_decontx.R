# Combine decontX results
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

# Load helper functions
source("./hubmap_single_cell/scRNA/scripts/single_cell_rna_helper_functions.R")

# Set directory paths
save_location <- "./hubmap_single_cell/scRNA_reproduce/initial_processing/"
individual_file_location <- "./hubmap_single_cell/scRNA_reproduce/initial_processing/new_decontx/"

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################

# Load and combine decontX results
setwd(individual_file_location)
decontXfiles <- list.files()[grepl("2.rds", list.files()) & grepl("decontX", list.files())]
project_names <- substr(decontXfiles, 8, nchar(decontXfiles)-5)
for (project_name in project_names){
	message(project_name)
	decontXresults <- readRDS(paste0("decontX", project_name, "2.rds"))
	new_counts <- decontXresults$decontXcounts
	colnames(new_counts) <- paste0(project_name, "_", colnames(new_counts))
	decontXcounts <- CreateSeuratObject(counts = new_counts)
	if (exists("full_counts")){
		full_counts <- merge(full_counts, y = decontXcounts)
	} else {
		full_counts <- decontXcounts
	}
}

setwd(save_location)
saveRDS(full_counts, "fulldecontxobj.rds")

