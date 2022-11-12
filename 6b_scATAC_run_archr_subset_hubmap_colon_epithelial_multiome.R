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
new_project_save_name = "HuBMAP_colon_epithelial_cells_multiome"
subscript = "epithelial_colon"

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 0) Load and subset previously defined archr project
if (0 %in% execute_steps){
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
        (proj$Location %in% c("Sigmoid","Transverse","Descending","Ascending")) & 
        (proj$CellTypeRNA %in% epi_cell_types) &
        (proj$Sample %ni% c("B009-A-001")))
    cellsSample <- proj$cellNames[idxSample]

    proj_stromal <- subsetArchRProject(
      ArchRProj = proj,
      cells = cellsSample,
      outputDirectory = new_project_save_name, dropCells = FALSE
    )

    saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)
} else {
    proj_stromal <- loadArchRProject(path = new_project_save_name)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 1) LSI Projection and Clustering
if (1 %in% execute_steps){
    library(parallel)
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
        groupBy = "Sample", force = TRUE
    )
    proj_stromal <- addUMAP(
        ArchRProj = proj_stromal, 
        reducedDims = paste("Harmony", subscript, sep = ""), 
        name = paste("UMAPHarmony", subscript, sep = ""), 
        nNeighbors = 50, 
        minDist = 0.3, 
        metric = "cosine", force = TRUE
    )
    proj_stromal <- addClusters(
        input = proj_stromal,
        reducedDims = paste("Harmony", subscript, sep = ""),
        method = "Seurat",
        name = paste("ClustersHarmony", subscript, sep = ""),
        resolution = 0.5, force=TRUE, nOutlier = 30, seed = 1, maxClusters = 40
    )
    p7 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = paste("ClustersHarmony", subscript, sep = ""), embedding = paste("UMAPHarmony", subscript, sep = ""))
    plotPDF(p7, name = paste(paste("Plot-UMAP2Harmony-Clusters", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)

    p3 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Sample", embedding = paste("UMAPHarmony", subscript, sep = ""))
    p5 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Location", embedding = paste("UMAPHarmony", subscript, sep = ""))
    p6 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = "Donor", embedding = paste("UMAPHarmony", subscript, sep = ""))
    p7 <- plotEmbedding(ArchRProj = proj_stromal, colorBy = "cellColData", name = paste("ClustersHarmony", subscript, sep = ""), embedding = paste("UMAPHarmony", subscript, sep = ""))
    plotPDF(p3,p5,p6,p7, name = paste(paste("Plot-UMAP2Harmony-Sample-Location-Donor-Harmony", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)
    saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)

    pal <- c()
    pal <- c("#371377")
    names(pal) <- "Best4+ Enterocytes"
    pal["CyclingTA 1"] <- "grey35"
    pal["CyclingTA 2"] <- "#cccccc"
    pal["D Cells"] <- "#FF0080"
    pal["Enterochromaffin"] <- "#FF0080"
    pal["Enterocytes"] <- "#9E0142"
    pal["EnteroendocrineUn"] <- "#FF0080"
    pal["EnteroendocrineUn 1"] <- "#FF0080"
    #pal["Exocrine"] <- "#88CCEE"
    #pal["Enteroendocrine2"] <- "#F59899"
    #pal["Enteroendocrine3"] <- "#E11E26"
    #pal["Epithelial"] <- "#238B45"
    pal["Goblet"] <- "#DC494C"
    pal["I Cells"] <- "#FF0080"
    pal["Immature Enterocytes"] <- "#7700FF"
    pal["Immature Goblet"] <- "#FAD510"
    #pal["K Cells"] <- "#FF0080"
    pal["L Cells"] <- "#FF0080"
    pal["Mo Cells"] <- "#FF0080"
    pal["NEUROG3high"] <- "#FF0080"
    #pal["Paneth"] <- "#F88D51"
    #pal["Secretory TA"] <- "#88CFA4"
    pal["S Cells"] <- "#FF0080"
    #pal["Secretory Specialized MUC5B+"] <- "#88CCEE"
    #pal["Secretory Specialized MUC6+"] <- "#88CCEE"
    #pal["Secretory Unknown 2"] <- "#2a7185"
    pal["Stem"] <- "#02401B"
    pal["TA1"] <- "#079290"
    pal["TA2"] <- "#A2A475"
    pal["Tuft"] <- "#0AD7D3"
    #pal["Unknown"] <- "#046C9A"
    p1 <- plotEmbedding(proj_stromal, colorBy = "cellColData", name = "CellTypeRNA", embedding = paste0("UMAPHarmony", subscript), pal = pal, size = 0.5)
    p2 <- plotEmbedding(proj_stromal, colorBy = "cellColData", name = "CellTypeRNA", embedding = paste0("UMAP", subscript), pal = pal)
    plotPDF(p1,p2, name = "Plot-UMAP-CellType.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)

    pal <- c()
    pal["CyclingTA 1"] <- "grey35"
    pal["CyclingTA 2"] <- "#cccccc"
    pal["Enterocytes"] <- "#9E0142"
    pal["Immature Enterocytes"] <- "#7700FF"
    pal["Stem"] <- "#02401B"
    pal["TA1"] <- "#079290"
    pal["TA2"] <- "#A2A475"

    p1 <- plotEmbedding(proj_stromal, colorBy = "cellColData", name = "CellTypeRNA", embedding = paste0("UMAPHarmony", subscript), pal = pal, size = 0.5)
    p2 <- plotEmbedding(proj_stromal, colorBy = "cellColData", name = "CellTypeRNA", embedding = paste0("UMAP", subscript), pal = pal)
    plotPDF(p1,p2, name = "Plot-UMAP-CellType_Abs_only.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)

}


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 2) Add RNA data to multiome project
colon <- readRDS(paste0(parent_directory, "scRNA/epithelial/colon/decontx_sctransform_cca/colon_clustered.rds"))

# add rna data
samples <- unique(c(
    paste0(colon$orig.ident)))

non_multiome <- c(samples[grepl("B001", samples)],
  samples[grepl("B004", samples)],
  samples[grepl("B005", samples)], "B001-A-302", "B004-A-004-R2")

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

rnaRowRanges <- readRDS("10x_RNA_row_ranges.rds")

colon.data.full <- colon.data.full[rownames(colon.data.full) %in% names(rnaRowRanges),]
library(stringr)
colnames(colon.data.full) <- str_replace(colnames(colon.data.full), "_", "#")

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
# 4) Add Peaks
if (4 %in% execute_steps){
    proj_ref <- loadArchRProject(path = "HuBMAP_epithelial_cells_multiome_final")
    peakSet <- getPeakSet(proj_ref)
    proj_stromal <- addPeakSet(ArchRProj = proj_stromal, peakSet = peakSet, force = TRUE)

    # Add Peak Matrix
    proj_stromal <- addPeakMatrix(ArchRProj = proj_stromal, force = TRUE)

    saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 6) Peak to gene linkages
proj_stromal <- addPeak2GeneLinks(
    ArchRProj = proj_stromal,
    reducedDims = paste0("Harmony", subscript),
    useMatrix = "GeneExpressionMatrix",
)
saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 7) Add deviations Matrix
motifPWMs <- readRDS("Vierstra-Human-Motifs.rds")
proj_stromal <- addMotifAnnotations(proj_stromal, motifPWMs = motifPWMs, name = "Vierstra", force = TRUE)
proj_stromal <- addBgdPeaks(proj_stromal)
proj_stromal <- addDeviationsMatrix(proj_stromal, peakAnnotation = "Vierstra", force = TRUE)
saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 9) Trajectory Analysis
trajectory <- c("Stem", "TA2", "TA1", "Immature Enterocytes", "Enterocytes")
trajectory <- c("C4", "C5", "C1")
proj_stromal <- addTrajectory(
    ArchRProj = proj_stromal, 
    name = "AbsorptiveU", 
    groupBy = paste("ClustersHarmony", subscript, sep = ""),#"CellTypeRNA",
    preFilterQuantile = 0.85,
    postFilterQuantile = 0.85,
    trajectory = trajectory, 
    reducedDims = paste0("Harmony", subscript),
    #embedding = paste0("UMAPHarmony", subscript), 
    force = TRUE
)

p <- plotTrajectory(proj_stromal, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmony", subscript), colorBy = "cellColData", name = "AbsorptiveU")
plotPDF(p, name = "Plot-AbsorptiveU-Traj-UMAPHarmony.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)

proj_stromal <- addImputeWeights(proj_stromal, reducedDims = paste0("Harmony", subscript))

p1 <- plotTrajectory(proj_stromal, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmony", subscript), colorBy = "GeneScoreMatrix", name = "ASCL2", continuousSet = "horizonExtra")
plotPDF(p1, name = "Plot-AbsorptiveU-Genes-Traj-UMAPHarmony.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 5, height = 5)

trajMM  <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU", useMatrix = "VierstraMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajPM  <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p2, p4, name = "Plot-AbsorptiveU-Traj-Heatmaps.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 6, height = 8)

trajPM  <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"), maxFeatures = 200000)
plotPDF(p4, name = "Plot-AbsorptiveU-Traj-Heatmaps_200K.pdf", ArchRProj = proj_stromal, addDOC = FALSE, width = 6, height = 8)

saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)


