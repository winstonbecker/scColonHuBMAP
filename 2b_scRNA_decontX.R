# run decontX for ambient RNA correction

.libPaths("./libraries/R_LIBS_4p1p2/")
library(celda)
library(stringr)
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
library(SingleCellExperiment)

###############################################################################################################################
###############################################################################################################################
# Load Args
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
j <- as.integer(args[2])

# Define variables
setwd("./hubmap_single_cell/scRNA_reproduce/initial_processing/")
scRNA_sets <- list(c("B001-A-001","B001-A-006","B001-A-101","B001-A-201", "B001-A-301","B001-A-401","B001-A-406", "B001-A-501"),
                    c("B004-A-004", "B004-A-008", "B004-A-104", "B004-A-204", "B004-A-304","B004-A-504","B004-A-404","B004-A-408"),#,"B004-A-404","B004-A-408"
                    c("B005-A-001","B005-A-002","B005-A-101","B005-A-201","B005-A-301","B005-A-401","B005-A-402","B005-A-501"), #,"B005-A-501"
                    c("B006-A-201", "B006-A-401", "B006-A-402","B006-A-201-R2", "B006-A-101", "B006-A-001", "B006-A-501", "B006-A-301"),
                    c("Bei_CO1-1_B010-A-101", "Bei_CO1-2_B010-A-201", "Bei_CO1-3_B010-A-405", "Bei_CO1-4_B010-A-501", "Bei_CO1-5_B010-A-001","Bei_CO1-6_B010-A-002", "Bei_CO1-7_B010-A-301", "Bei_CO1-8_B010-A-401"), 
                    c("Bei_CO2-1_B011-A-001", "Bei_CO2-10_B012-A-002","Bei_CO2-11_B012-A-301", "Bei_CO2-12_B012-A-401", "Bei_CO2-13_B012-A-101", "Bei_CO2-14_B012-A-201","Bei_CO2-17_B008-A-002", "Bei_CO2-18_B008-A-101"),
                    c("Bei_CO2-19_B008-A-201", "Bei_CO2-2_B011-A-002", "Bei_CO2-20_B008-A-401", "Bei_CO2-21_B008-A-402", "Bei_CO2-22_B008-A-501", "Bei_CO2-23_B012-A-405", "Bei_CO2-24_B012-A-501", "Bei_CO2-25_B009-A-001"), 
                    c("Bei_CO2-26_B008-A-001", "Bei_CO2-27_B009-A-301", "Bei_CO2-28_B009-A-401", "Bei_CO2-29_B009-A-101", "Bei_CO2-3_B011-A-301", "Bei_CO2-30_B008-A-301", "Bei_CO2-31_B009-A-405"), #, "Bei_CO2-32_B009-A-501" 
                    c("B004-A-404-R2", "B006-A-002", "Bei_CO2-4_B011-A-401", "Bei_CO2-5_B011-A-101","Bei_CO2-6_B011-A-201", "Bei_CO2-7_B011-A-405", "Bei_CO2-8_B011-A-501", "Bei_CO2-9_B012-A-001"),
                    c("B004-A-104-R2", "B004-A-408-R2", "B004-A-504-R2", "B005-A-501-R2"))

# Define functions
run_decontx <- function(project_name){
  # Function for seurat object creation and basic qc and filtering and running of doublet finder
  # project_name: name that will be used for project in the seurat object and that will be used when saving plots
  project_name <- rev(str_split(project_name, "_")[[1]])[1]


  filtered_counts <- as.matrix(readRDS(paste0("filtered_counts_matrix", project_name, ".rds")))
  storage.mode(filtered_counts) <- "integer"
  #filtered_counts <- filtered_counts[rowSums(filtered_counts)>0,]

  decontXresults <- decontX(filtered_counts)#, background = raw_counts)

  saveRDS(decontXresults, paste0("decontX", project_name, "2.rds"))

  # sce <- SingleCellExperiment(list(counts = filtered_counts))
  # sce <- decontX(sce)

  # umap <- reducedDim(sce, "decontX_UMAP")
  # pdf(paste0("./", project_name, "_decontX", "_UMAP.pdf"))
  # print(plotDimReduceCluster(x = sce$decontX_clusters,
  #     dim1 = umap[, 1], dim2 = umap[, 2]))
  # dev.off()

  # pdf(paste0("./", project_name, "_decontX", "_contamination_UMAP.pdf"))
  # print(plotDecontXContamination(sce))
  # dev.off()
}

###############################################################################################################################
###############################################################################################################################
# Run
samples <- scRNA_sets[[j]]
run_decontx(samples[i])


