# Set things up
.libPaths("./libraries/R_LIBS_4p1p2/")

# Load packages
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
`%notin%` <- Negate(`%in%`)

#Set/Create Working Directory to Folder
parent_directory <- "/hubmap_single_cell/"
setwd(paste0(parent_directory, "scATAC/projects/"))

#Set Threads to be used
addArchRThreads()

# Things to set for subseting projects. Also find and replace proj name for subsample with desired project name (or similar)
new_project_save_name = "HuBMAP_epithelial_cells_multiome_final"
subscript = "epithelial"

set.seed(1)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 0) Load and subset previously defined archr project

proj <- loadArchRProject("all_hubmap_cells_filtered")

cell_types_scrna <- read.table("scrna_cell_types.tsv")
library(stringr)
rownames(cell_types_scrna) <- str_replace(rownames(cell_types_scrna), "_", "#")
multiome <- c(rownames(cell_types_scrna)[grepl("B006", rownames(cell_types_scrna))],
  rownames(cell_types_scrna)[grepl("B008", rownames(cell_types_scrna))],
  rownames(cell_types_scrna)[grepl("B009", rownames(cell_types_scrna))],
  rownames(cell_types_scrna)[grepl("B010", rownames(cell_types_scrna))],
  rownames(cell_types_scrna)[grepl("B011", rownames(cell_types_scrna))],
  rownames(cell_types_scrna)[grepl("B012", rownames(cell_types_scrna))])
cell_types_scrna <- cell_types_scrna[rownames(cell_types_scrna) %in% multiome,, drop = FALSE]
cell_types_scrna <- cell_types_scrna[rownames(cell_types_scrna) %in% rownames(getCellColData(proj)),, drop = FALSE]

proj <- addCellColData(ArchRProj = proj, data = paste0(cell_types_scrna$CellType), cells = paste0(rownames(cell_types_scrna)), name = "CellTypeRNA", force = TRUE)

epi_cell_types <- c("Immature Enterocytes",
                "Enterocytes",
                "CyclingTA 1",
                "TA2",
                "Epithelial",
                "Goblet",
                "Best4+ Enterocytes",
                "Immature Goblet",
                "TA1",
                "Stem",
                "Enteroendocrine",
                "Paneth",
                "Tuft",
                "CyclingTA 2",
                "Secretory Specialized MUC6+",
                "Unknown",
                "EnteroendocrineUn 1",
                "Enterochromaffin",
                "NEUROG3high",
                "L Cells",
                "I Cells",
                "EnteroendocrineUn",
                "Mo Cells",
                "D Cells",
                "K Cells",
                "S Cells",
                "Exocrine",
                "Secretory Specialized MUC5B+")

# Can start with all cells filtered project and then subset
# Define multiome subset to explore
idxSample <- BiocGenerics::which((proj$Sample %in% unique(getCellColData(proj)$Sample)[order(unique(getCellColData(proj)$Sample))][27:71]) & 
    (proj$CellTypeRNA %in% epi_cell_types) & 
    (proj$Sample %ni% c("B009-A-001", "B009-A-301")))

cellsSample <- proj$cellNames[idxSample]

proj_stromal <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsSample,
  outputDirectory = new_project_save_name, dropCells = FALSE
)

saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# LSI Projection and Clustering
proj_stromal <- addIterativeLSI(
    ArchRProj = proj_stromal,
    useMatrix = "TileMatrix", 
    name = paste("IterativeLSI", subscript, sep = ""), 
    iterations = 3, 
    clusterParams = list(
        resolution = c(0.1, 0.2), 
        sampleCells = NULL, 
        n.start = 10
    ), 
    varFeatures = 20000, sampleCellsPre = NULL,
    dimsToUse = 1:30, force = TRUE
)
proj_stromal <- addUMAP(
    ArchRProj = proj_stromal, 
    reducedDims = paste("IterativeLSI", subscript, sep = ""), 
    name = paste("UMAP", subscript, sep = ""), 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", force=TRUE
)

p1 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Sample", embedding = paste("UMAP", subscript, sep = ""))
p3 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Donor", embedding = paste("UMAP", subscript, sep = ""))
p4 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Location", embedding = paste("UMAP", subscript, sep = ""))
plotPDF(p1,p3,p4, name = paste(paste("Plot-UMAP-Sample-Clusters-Donor-Location", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)

proj_stromal <- addHarmony(
    ArchRProj = proj_stromal,
    reducedDims = paste("IterativeLSI", subscript, sep = ""),
    name = paste("Harmony", subscript, sep = ""),
    groupBy = "Donor", force = TRUE
)
proj_stromal <- addUMAP(
    ArchRProj = proj_stromal, 
    reducedDims = paste("Harmony", subscript, sep = ""), 
    name = paste("UMAPHarmony", subscript, sep = ""), 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", force = TRUE
)
proj_stromal <- addClusters(
    input = proj_stromal,
    reducedDims = paste("Harmony", subscript, sep = ""),
    method = "Seurat",
    name = paste("ClustersHarmony", subscript, sep = ""),
    resolution = 1.2, force=TRUE, nOutlier = 50, seed = 1
)

p3 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Sample", embedding = paste("UMAPHarmony", subscript, sep = ""))
p5 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Location", embedding = paste("UMAPHarmony", subscript, sep = ""))
p6 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Donor", embedding = paste("UMAPHarmony", subscript, sep = ""))
plotPDF(p3,p5,p6, name = paste(paste("Plot-UMAP2Harmony-Sample-Location-Donor-Harmony", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)

p7 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "CellTypeRNA", embedding = paste("UMAPHarmony", subscript, sep = ""))
plotPDF(p7, name = paste(paste("Plot-UMAP2Harmony-Sample-Location-Donor-Harmony", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)

proj_stromal <- addHarmony(
    ArchRProj = proj_stromal,
    reducedDims = paste("IterativeLSI", subscript, sep = ""),
    name = paste("HarmonySample", subscript, sep = ""),
    groupBy = "Sample", force = TRUE
)
proj_stromal <- addUMAP(
    ArchRProj = proj_stromal, 
    reducedDims = paste("HarmonySample", subscript, sep = ""), 
    name = paste("UMAPHarmonySample", subscript, sep = ""), 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", force = TRUE
)

p3 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Sample", embedding = paste("UMAPHarmonySample", subscript, sep = ""))
p5 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Location", embedding = paste("UMAPHarmonySample", subscript, sep = ""))
p6 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Donor", embedding = paste("UMAPHarmonySample", subscript, sep = ""))
#p7 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "CellType", embedding = paste("UMAPHarmonySample", subscript, sep = ""))
plotPDF(p3,p5,p6,p7, name = paste(paste("Plot-UMAPHarmonySample-Sample-Location-Donor-Harmony", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Add RNA data to multiome project
colon <- readRDS(paste0(parent_directory, "scRNA/epithelial/colon/decontx_sctransform_cca/colon_clustered.rds"))
ileum <- readRDS(paste0(parent_directory, "scRNA/epithelial/ileum/decontx_sctransform_cca/ileum.rds"))
jejunum <- readRDS(paste0(parent_directory, "scRNA/epithelial/jejunum/decontx_sctransform_cca/jejunum_clustered.rds"))
duodenum <- readRDS(paste0(parent_directory, "scRNA/epithelial/duodenum/decontx_sctransform_cca/duodenum_clustered.rds"))
enteroendocrine <- readRDS(paste0(parent_directory, "scRNA/epithelial/enteroendocrine/decontx_sctransform_cca/clustered_annotated_enteroendocrine_cells.rds"))
secretory_special <- readRDS(paste0(parent_directory, "scRNA/epithelial/specialized_secretory/decontx_sctransform_cca/specialized_secretory_clustered.rds"))

DefaultAssay(colon) <- "decontXcounts"
colon <- DietSeurat(colon, assays = c("RNA", "decontXcounts"))

DefaultAssay(ileum) <- "decontXcounts"
ileum <- DietSeurat(ileum, assays = c("RNA", "decontXcounts"))

DefaultAssay(jejunum) <- "decontXcounts"
jejunum <- DietSeurat(jejunum, assays = c("RNA", "decontXcounts"))

DefaultAssay(duodenum) <- "decontXcounts"
duodenum <- DietSeurat(duodenum, assays = c("RNA", "decontXcounts"))

DefaultAssay(enteroendocrine) <- "decontXcounts"
enteroendocrine <- DietSeurat(enteroendocrine, assays = c("RNA", "decontXcounts"))

DefaultAssay(secretory_special) <- "decontXcounts"
secretory_special <- DietSeurat(secretory_special, assays = c("RNA", "decontXcounts"))

colon <- subset(x = colon, subset = CellType %ni% c("Enteroendocrine", "Secretory Specialized MUC6+"))
ileum <- subset(x = ileum, subset = CellType %ni% c("Enteroendocrine", "Secretory Specialized MUC6+"))
jejunum <- subset(x = jejunum, subset = CellType %ni% c("Enteroendocrine", "Secretory Specialized MUC6+"))
duodenum <- subset(x = duodenum, subset = CellType %ni% c("Enteroendocrine", "Secretory Specialized MUC6+"))
colon <- merge(colon, y = c(duodenum, ileum, jejunum, enteroendocrine, secretory_special))

# add rna data
samples <- unique(c(
    paste0(colon$orig.ident)))

multiome <- c(samples[grepl("B006", samples)],
  samples[grepl("B008", samples)],
  samples[grepl("B009", samples)],
  samples[grepl("B010", samples)],
  samples[grepl("B011", samples)],
  samples[grepl("B012", samples)])

DefaultAssay(colon) <- "decontXcounts"
colon_diet <- DietSeurat(colon, assays = c("decontXcounts"))
colon_diet <- subset(colon_diet, subset = orig.ident %in% multiome)
colon.data.full <- GetAssayData(object = colon_diet, slot = "counts")

seRNA <- import10xFeatureMatrix(
    input = paste0(c("./data/B008-A-402/outs/raw_feature_bc_matrix.h5")),
    names = c("B008-A-402")
)
rnaRowRanges <- rowRanges(seRNA)
saveRDS(rnaRowRanges, "10x_RNA_row_ranges.rds")

colon.data.full <- colon.data.full[rownames(colon.data.full) %in% names(rnaRowRanges),]
library(stringr)
colnames(colon.data.full) <- str_replace(colnames(colon.data.full), "_", "#")
#
# colon.data.full.init <- GetAssayData(object = colon_diet, slot = "counts")
# colon.data.full.init <- colon.data.full.init[rownames(colon.data.full.init) %ni% names(rnaRowRanges),]

rowRangesSE <- rnaRowRanges[names(rnaRowRanges) %in% rownames(colon.data.full),]
colon.data.full <- colon.data.full[names(rowRangesSE),] # make sure they have the same order

seRNA<-SummarizedExperiment(assays=list(counts=colon.data.full), rowRanges= rowRangesSE) # create summarized experiment for adding to arrow
library(parallel)
proj_stromal <- addGeneExpressionMatrix(input = proj_stromal, seRNA = seRNA, force = TRUE)

# save project
saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Call Peaks
cellTypes <- getCellColData(proj_stromal)[,c("CellTypeRNA"), drop = FALSE]
cellTypes[cellTypes$CellTypeRNA == "K Cells",] <- "Enteroendocrine"
cellTypes[cellTypes$CellTypeRNA == "L Cells",] <- "Enteroendocrine"
cellTypes[cellTypes$CellTypeRNA == "Mo Cells",] <- "Enteroendocrine"
cellTypes[cellTypes$CellTypeRNA == "D Cells",] <- "Enteroendocrine"
cellTypes[cellTypes$CellTypeRNA == "I Cells",] <- "Enteroendocrine"
cellTypes[cellTypes$CellTypeRNA == "S Cells",] <- "Enteroendocrine"
cellTypes[cellTypes$CellTypeRNA == "NEUROG3high",] <- "Enteroendocrine"
cellTypes[cellTypes$CellTypeRNA == "Enterochromaffin",] <- "Enteroendocrine"
cellTypes[cellTypes$CellTypeRNA == "EnteroendocrineUn 1",] <- "Enteroendocrine"
cellTypes[cellTypes$CellTypeRNA == "EnteroendocrineUn",] <- "Enteroendocrine"
cellTypes[cellTypes$CellTypeRNA == "Secretory Specialized MUC5B+",] <- "Secretory Specialized MUC6+"
cellTypes[cellTypes$CellTypeRNA == "Exocrine",] <- "Secretory Specialized MUC6+"
cellTypes[cellTypes$CellTypeRNA == "Unknown",] <- "Secretory Specialized MUC6+"
cellLocations <- getCellColData(proj_stromal)[,c("Location"), drop = FALSE]
cellLocations[cellLocations$Location == "Transverse",] <- "Colon"
cellLocations[cellLocations$Location == "Sigmoid",] <- "Colon"
cellLocations[cellLocations$Location == "Descending",] <- "Colon"
cellLocations[cellLocations$Location == "Ascending",] <- "Colon"
cellLocations[cellLocations$Location == "Mid-jejunum",] <- "Proximal-jejunum"
cellTypes$CellTypeRNA <- paste0(cellLocations$Location, "-", cellTypes$CellTypeRNA)
cellTypes[cellTypes$CellTypeRNA == "Duodenum-Tuft",] <- "Proximal-jejunum-Tuft"
cellTypes[cellTypes$CellTypeRNA == "Ileum-Tuft",] <- "Proximal-jejunum-Tuft"
proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste0(cellTypes$CellTypeRNA), cells = rownames(cellTypes), name = "CellTypeLocation", force = TRUE)

set.seed(1)
proj_stromal <- addGroupCoverages(ArchRProj = proj_stromal, groupBy = "CellTypeLocation", force = TRUE)
pathToMacs2 <- findMacs2()

#Call Reproducible Peaks w/ Macs2
proj_stromal <- addReproduciblePeakSet(
    ArchRProj = proj_stromal, groupBy = "CellTypeLocation", force = TRUE, 
    pathToMacs2 = pathToMacs2, maxPeaks = 250000
)

#Add Peak Matrix
proj_stromal <- addPeakMatrix(ArchRProj = proj_stromal, force = TRUE)

saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Add deviations matrix
system("wget https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra-Human-Motifs.rds")
motifPWMs <- readRDS("Vierstra-Human-Motifs.rds")
proj_stromal <- addMotifAnnotations(proj_stromal, motifPWMs = motifPWMs, name = "Vierstra")
proj_stromal <- addBgdPeaks(proj_stromal)
proj_stromal <- addDeviationsMatrix(proj_stromal, peakAnnotation = "Vierstra", force = TRUE)
saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)



