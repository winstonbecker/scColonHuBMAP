# Script for initial qc and arrow creation for scATAC data from HuBMAP normal intestine project. 

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Load packages
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
`%notin%` <- Negate(`%in%`)

#Load Genome Annotations
addArchRGenome("hg38")

#Set Threads to be used
addArchRThreads()

#Set/Create Working Directory to Folder
setwd("./hubmap_single_cell/scATAC/projects/")

singleome_fragments_files_directory <- "./scHuBMAP_HTAN/scATAC/fragments_files/"
multiome_cell_ranger_outputs_directory <- "./scHuBMAP_HTAN/hubmap_multiome_plus_redos/cell_ranger_outputs/"

################################################################################################################
#..............................................................................................................#
################################################################################################################
# Make arrow files
# Define path to the individaul fragments files

inputFiles1K <- c(
  "B006-A-001" = paste0(multiome_cell_ranger_outputs_directory, "B006-A-001/outs/B006-A-001_atac_fragments.tsv.gz"),
  "B006-A-101" = paste0(multiome_cell_ranger_outputs_directory, "B006-A-101/outs/B006-A-101_atac_fragments.tsv.gz"),
  "B006-A-201-R2" = paste0(multiome_cell_ranger_outputs_directory, "B006-A-201-R2/outs/B006-A-201-R2_atac_fragments.tsv.gz"),
  "B006-A-301" = paste0(multiome_cell_ranger_outputs_directory, "B006-A-301/outs/B006-A-301_atac_fragments.tsv.gz"),
  "B006-A-501" = paste0(multiome_cell_ranger_outputs_directory, "B006-A-501/outs/B006-A-501_atac_fragments.tsv.gz")
)
inputFiles1p5K <- c(
  "B001-A-302" = paste0(singleome_fragments_files_directory, "HuBMAP/IN_HTAN_PAPER/B001-A-302-D_fragments.tsv.gz"),
  "B001-A-401" = paste0(singleome_fragments_files_directory, "HuBMAP/IN_HTAN_PAPER/B001-A-401-D_fragments.tsv.gz"),
  "B001-A-406" = paste0(singleome_fragments_files_directory, "HuBMAP/IN_HTAN_PAPER/B001-A-406-D_fragments.tsv.gz"),
  "B001-A-501" = paste0(singleome_fragments_files_directory, "HuBMAP/IN_HTAN_PAPER/B001-A-501-D_fragments.tsv.gz")
)
inputFiles2K <- c(
  "B001-A-001" = paste0(singleome_fragments_files_directory, "HuBMAP/B001-A-001-D_20200715_fragments.tsv.gz"),
  "B001-A-006" = paste0(singleome_fragments_files_directory, "HuBMAP/B001-A-006-D_20200715_fragments.tsv.gz"),
  "B001-A-101" = paste0(singleome_fragments_files_directory, "HuBMAP/B001-A-101-D_20200804_fragments.tsv.gz"),
  "B004-A-004" = paste0(singleome_fragments_files_directory, "HuBMAP/IN_HTAN_PAPER/B004-A-004-D_20200715_fragments.tsv.gz"), 
  "B004-A-404" = paste0(singleome_fragments_files_directory, "HuBMAP/B004-A-404-D_20200715_fragments.tsv.gz"), 
  "B008-A-001" = paste0(multiome_cell_ranger_outputs_directory, "B008-A-001/outs/atac_fragments.tsv.gz")
)
inputFiles3K <- c(
  "B001-A-201" = paste0(singleome_fragments_files_directory, "HuBMAP/B001-A-201-D_08042020_fragments.tsv.gz"),
  "B001-A-301" = paste0(singleome_fragments_files_directory, "HuBMAP/IN_HTAN_PAPER/B001-A-301-D_20200804_fragments.tsv.gz"),
  "B004-A-404-R2" = paste0(singleome_fragments_files_directory, "HuBMAP/B004-A-404-D_20200817_fragments.tsv.gz"),
  "B004-A-408" = paste0(singleome_fragments_files_directory, "HuBMAP/B004-A-408-D_20200715_fragments.tsv.gz"),
  "B004-A-504" = paste0(singleome_fragments_files_directory, "HuBMAP/B004-A-504-D_20200715_fragments.tsv.gz"),
  "B005-A-301" = paste0(singleome_fragments_files_directory, "HuBMAP/B005-A-301-D_20200917_fragments.tsv.gz"),
  "B005-A-501" = paste0(singleome_fragments_files_directory, "HuBMAP/B005-A-501-D_20200917_fragments.tsv.gz"),
  "B005-A-001" = paste0(singleome_fragments_files_directory, "HuBMAP/B005-A-001_fragments.tsv.gz"),
  "B005-A-002" = paste0(singleome_fragments_files_directory, "HuBMAP/B005-A-002_fragments.tsv.gz"),
  "B005-A-101" = paste0(singleome_fragments_files_directory, "HuBMAP/B005-A-101_fragments.tsv.gz"),
  "B005-A-201" = paste0(singleome_fragments_files_directory, "HuBMAP/B005-A-201_fragments.tsv.gz"),
  "B005-A-401" = paste0(singleome_fragments_files_directory, "HuBMAP/B005-A-401_fragments.tsv.gz"),
  "B005-A-402" = paste0(singleome_fragments_files_directory, "HuBMAP/B005-A-402_fragments.tsv.gz"),
  "B006-A-002" = paste0(multiome_cell_ranger_outputs_directory, "B006-A-002/outs/B006-A-002_atac_fragments.tsv.gz"),
  "B006-A-201" = paste0(multiome_cell_ranger_outputs_directory, "B006-A-201/outs/B006-A-201_atac_fragments.tsv.gz"),
  "B006-A-401" = paste0(multiome_cell_ranger_outputs_directory, "B006-A-401/outs/B006-A-401_atac_fragments.tsv.gz"),
  "B006-A-402" = paste0(multiome_cell_ranger_outputs_directory, "B006-A-402/outs/B006-A-402_atac_fragments.tsv.gz"),
  "B008-A-002" = paste0(multiome_cell_ranger_outputs_directory, "B008-A-002/outs/atac_fragments.tsv.gz"),
  "B008-A-101" = paste0(multiome_cell_ranger_outputs_directory, "B008-A-101/outs/atac_fragments.tsv.gz"),
  "B008-A-201" = paste0(multiome_cell_ranger_outputs_directory, "B008-A-201/outs/atac_fragments.tsv.gz"), 
  "B008-A-501" = paste0(multiome_cell_ranger_outputs_directory, "B008-A-501/outs/atac_fragments.tsv.gz"), 
  "B009-A-301" = paste0(multiome_cell_ranger_outputs_directory, "B009-A-301/outs/atac_fragments.tsv.gz"), 
  "B009-A-405" = paste0(multiome_cell_ranger_outputs_directory, "B009-A-405/outs/atac_fragments.tsv.gz"), 
  "B010-A-001" = paste0(multiome_cell_ranger_outputs_directory, "B010-A-001/outs/atac_fragments.tsv.gz"), 
  "B010-A-002" = paste0(multiome_cell_ranger_outputs_directory, "B010-A-002/outs/atac_fragments.tsv.gz"), 
  "B010-A-101" = paste0(multiome_cell_ranger_outputs_directory, "B010-A-101/outs/atac_fragments.tsv.gz"), 
  "B010-A-501" = paste0(multiome_cell_ranger_outputs_directory, "B010-A-501/outs/atac_fragments.tsv.gz"), 
  "B011-A-001" = paste0(multiome_cell_ranger_outputs_directory, "B011-A-001/outs/atac_fragments.tsv.gz"), 
  "B011-A-002" = paste0(multiome_cell_ranger_outputs_directory, "B011-A-002/outs/atac_fragments.tsv.gz"), 
  "B011-A-201" = paste0(multiome_cell_ranger_outputs_directory, "B011-A-201/outs/atac_fragments.tsv.gz"), 
  "B011-A-401" = paste0(multiome_cell_ranger_outputs_directory, "B011-A-401/outs/atac_fragments.tsv.gz"), 
  "B012-A-001" = paste0(multiome_cell_ranger_outputs_directory, "B012-A-001/outs/atac_fragments.tsv.gz"), 
  "B012-A-002" = paste0(multiome_cell_ranger_outputs_directory, "B012-A-002/outs/atac_fragments.tsv.gz"), 
  "B012-A-101" = paste0(multiome_cell_ranger_outputs_directory, "B012-A-101/outs/atac_fragments.tsv.gz"), 
  "B012-A-401" = paste0(multiome_cell_ranger_outputs_directory, "B012-A-401/outs/atac_fragments.tsv.gz"),
  "B009-A-501" = paste0(multiome_cell_ranger_outputs_directory, "B009-A-501/outs/atac_fragments.tsv.gz")
)
inputFiles4K <- c(
  "B004-A-004-R2" = paste0(singleome_fragments_files_directory, "HuBMAP/IN_HTAN_PAPER/B004-A-004-D_20200817_fragments.tsv.gz"), 
  "B004-A-204" = paste0(singleome_fragments_files_directory, "HuBMAP/IN_HTAN_PAPER/B004-A-204-D_20200702_fragments.tsv.gz"),
  "B004-A-304" = paste0(singleome_fragments_files_directory, "HuBMAP/B004-A-304-D_20200702_fragments.tsv.gz"),
  "B008-A-301" = paste0(multiome_cell_ranger_outputs_directory, "B008-A-301/outs/atac_fragments.tsv.gz"), 
  "B009-A-101" = paste0(multiome_cell_ranger_outputs_directory, "B009-A-101/outs/atac_fragments.tsv.gz"),
  "B010-A-201" = paste0(multiome_cell_ranger_outputs_directory, "B010-A-201/outs/atac_fragments.tsv.gz"), 
  "B010-A-401" = paste0(multiome_cell_ranger_outputs_directory, "B010-A-401/outs/atac_fragments.tsv.gz"), 
  "B011-A-101" = paste0(multiome_cell_ranger_outputs_directory, "B011-A-101/outs/atac_fragments.tsv.gz"), 
  "B011-A-301" = paste0(multiome_cell_ranger_outputs_directory, "B011-A-301/outs/atac_fragments.tsv.gz"), 
  "B011-A-501" = paste0(multiome_cell_ranger_outputs_directory, "B011-A-501/outs/atac_fragments.tsv.gz"), 
  "B012-A-501" = paste0(multiome_cell_ranger_outputs_directory, "B012-A-501/outs/atac_fragments.tsv.gz")
)
inputFiles5K <- c(
  "B004-A-008" = paste0(singleome_fragments_files_directory, "HuBMAP/IN_HTAN_PAPER/B004-A-008-D_20200817_fragments.tsv.gz"),
  "B010-A-301" = paste0(multiome_cell_ranger_outputs_directory, "B010-A-301/outs/atac_fragments.tsv.gz"), 
  "B011-A-405" = paste0(multiome_cell_ranger_outputs_directory, "B011-A-405/outs/atac_fragments.tsv.gz"), 
  "B012-A-405" = paste0(multiome_cell_ranger_outputs_directory, "B012-A-405/outs/atac_fragments.tsv.gz"), 
  "B008-A-402" = paste0(multiome_cell_ranger_outputs_directory, "B008-A-402/outs/atac_fragments.tsv.gz"), 
  "B009-A-001" = paste0(multiome_cell_ranger_outputs_directory, "B009-A-001/outs/atac_fragments.tsv.gz")
)
inputFiles6K <- c(
  "B012-A-201" = paste0(multiome_cell_ranger_outputs_directory, "B012-A-201/outs/atac_fragments.tsv.gz"), 
  "B008-A-401" = paste0(multiome_cell_ranger_outputs_directory, "B008-A-401/outs/atac_fragments.tsv.gz")
)
inputFiles10K <- c(
  "B010-A-405" = paste0(multiome_cell_ranger_outputs_directory, "B010-A-405/outs/atac_fragments.tsv.gz"),
  "B012-A-301" = paste0(multiome_cell_ranger_outputs_directory, "B012-A-301/outs/atac_fragments.tsv.gz")
)
allInput <- c(inputFiles1K, inputFiles1p5K, inputFiles2K, inputFiles3K, inputFiles4K, inputFiles5K, inputFiles6K, inputFiles10K)


ArrowFiles1K <- createArrowFiles(
  inputFiles = inputFiles1K,
  sampleNames = names(inputFiles1K), minFrags = 1000, minTSS = 5,
)
ArrowFiles1p5K <- createArrowFiles(
  inputFiles = inputFiles1p5K,
  sampleNames = names(inputFiles1p5K), minFrags = 1500, minTSS = 5,
)
ArrowFiles2K <- createArrowFiles(
  inputFiles = inputFiles2K,
  sampleNames = names(inputFiles2K), minFrags = 2000, minTSS = 5,
)
ArrowFiles3K <- createArrowFiles(
  inputFiles = inputFiles3K,
  sampleNames = names(inputFiles3K), minFrags = 3000, minTSS = 5,
)
ArrowFiles4K <- createArrowFiles(
  inputFiles = inputFiles4K,
  sampleNames = names(inputFiles4K), minFrags = 4000, minTSS = 5,
)
ArrowFiles5K <- createArrowFiles(
  inputFiles = inputFiles5K,
  sampleNames = names(inputFiles5K), minFrags = 5000, minTSS = 5,
)
ArrowFiles6K <- createArrowFiles(
  inputFiles = inputFiles6K,
  sampleNames = names(inputFiles6K), minFrags = 6000, minTSS = 5,
)
ArrowFiles10K <- createArrowFiles(
  inputFiles = inputFiles10K,
  sampleNames = names(inputFiles10K), minFrags = 10000, minTSS = 5,
)

ArrowFiles1K <- paste0(names(inputFiles1K), ".arrow")
ArrowFiles1p5K <- paste0(names(inputFiles1p5K), ".arrow")
ArrowFiles2K <- paste0(names(inputFiles2K), ".arrow")
ArrowFiles3K <- paste0(names(inputFiles3K), ".arrow")
ArrowFiles4K <- paste0(names(inputFiles4K), ".arrow")
ArrowFiles5K <- paste0(names(inputFiles5K), ".arrow")
ArrowFiles6K <- paste0(names(inputFiles6K), ".arrow")
ArrowFiles10K <- paste0(names(inputFiles10K), ".arrow")

doubScores <- addDoubletScores(ArrowFiles1K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
doubScores <- addDoubletScores(ArrowFiles1p5K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
doubScores <- addDoubletScores(ArrowFiles2K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
doubScores <- addDoubletScores(ArrowFiles3K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
doubScores <- addDoubletScores(ArrowFiles4K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
doubScores <- addDoubletScores(ArrowFiles5K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
doubScores <- addDoubletScores(ArrowFiles6K, k = 10, knnMethod = "UMAP", LSIMethod = 1)
doubScores <- addDoubletScores(ArrowFiles10K, k = 10, knnMethod = "UMAP", LSIMethod = 1)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Create ArchR project
set.seed(1)
ArrowFiles1 <- c(ArrowFiles1K, ArrowFiles1p5K, ArrowFiles2K, ArrowFiles3K, ArrowFiles4K, ArrowFiles5K, ArrowFiles6K, ArrowFiles10K)
proj <- ArchRProject(
    ArrowFiles = c(ArrowFiles1), 
    outputDirectory = "all_hubmap_cells"
)
saveArchRProject(ArchRProj = proj, outputDirectory = "all_hubmap_cells", load = FALSE, overwrite = FALSE)
proj <- loadArchRProject("all_hubmap_cells")
proj <- filterDoublets(proj, filterRatio = 1.2)
write.table(rownames(getCellColData(proj)), "initial_post_filter_atac_cells.txt")

