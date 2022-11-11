# Script to analyze epithelial scRNA cells from HuBMAP intestine data
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
analysis_parent_folder <- paste0(parent_folder, "scRNA_reproduce/all_cells/decontx_norm_scale_harmony/") # folder for this analysis
setwd(analysis_parent_folder)

# load helper functions
source(paste0(parent_folder, "scRNA_reproduce/scripts/single_cell_rna_helper_functions.R"))

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# load and merge seurat objects containing all epithelial cells
colon <- readRDS(paste0(parent_folder, "scRNA/all_cells/decontx_norm_scale_harmony/colon_epithelial_all_samples_initial.rds"))
extra_cells <- readRDS(paste0(parent_folder, "scRNA_reproduce/stromal/decontx_norm_scale_harmony/stromal_inital_move_to_epithelial.rds"))
extra_cells1 <- readRDS(paste0(parent_folder, "scRNA_reproduce/stromal/decontx_norm_scale_harmony/stromal_inital_move_to_epithelial_R2.rds"))
extra_cells2 <- readRDS(paste0(parent_folder, "scRNA_reproduce/stromal/decontx_norm_scale_harmony/stromal_inital_move_to_epithelial_R3.rds"))
extra_cells3 <- readRDS(paste0(parent_folder, "scRNA_reproduce/immune/decontx_norm_scale_harmony/immune_inital_move_to_epithelail.rds"))
epithelial <- merge(colon, y = c(extra_cells, extra_cells1, extra_cells2, extra_cells3), project = "full_epithelial")
epithelial <- DietSeurat(epithelial)

# Subset for each compartment
duodenum <- DietSeurat(subset(epithelial, subset = Location == "Duodenum"))
saveRDS(duodenum, file = "duodenum_epithelial_samples_initial.rds")

jejunum <- DietSeurat(subset(epithelial, subset = Location %in% c("Proximal-jejunum", "Mid-jejunum")))
saveRDS(jejunum, file = "jejunum_epithelial_samples_initial.rds")

ileum <- DietSeurat(subset(epithelial, subset = Location == "Ileum"))
saveRDS(ileum, file = "ileum_epithelial_samples_initial.rds")

colon <- DietSeurat(subset(epithelial, subset = Location %in% c("Ascending", "Transverse", "Descending", "Sigmoid")))
saveRDS(colon, file = "colon_epithelial_samples_initial.rds")


