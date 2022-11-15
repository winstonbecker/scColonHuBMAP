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

# Things to set for subseting projects--find and replace proj name for subsample with desired project name (or similar)
subscript = "immune"
new_project_save_name <- "HuBMAP_immune_cells_multiome"
use_Regev_RNA <- FALSE

markerGenes  <- c(
    "PAX5", "MS4A1", "CD19", "IGLL5", "VPREB3", #B-Cell Trajectory
    "TPSAB1", "HDC", "CTSG", "CMA1", "KRT1", "IL1RAPL1", "GATA2", #Mast
    "SERPINA9", "HRK", "HTR3A", "TCL6", "CD180", "FCRLA", #GC
    "CMA1", "IL1RAPL1", "CD69", #CD69+ Mast
    "KRT1", #CD69- Mast
    "CD207", #DC2
    "KLRF1", "SH2D1B", "SH2D1B", #NKs
    "SSR4", "IGLL5", "IGLL1", "AMPD1",#Plasma
    "CD14", "CLEC9A", "FCGR1A", "LILRB2", "CD209", "CD1E", #Monocytes
    "S100A8", "S100A9", # Inflammatory Monocytes
    "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "TBX21", "IL7R", "CD4", "CD2", #TCells
    "BATF","TNFRSF4", "FOXP3","CTLA4","LAIR2" # Tregs
    )
# Define larger dataframe of known RNA markers
RNA_markers <- read.csv(file = '/oak/stanford/groups/wjg/wbecker/other/scATAC/HuBMAP_HTAN_ENCODE_Only/RegevLabMarkers/S2_Regev_Cell_immune.csv')

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 0) Load and subset previously defined archr project
if (0 %in% execute_steps){
    # Load previously defined archr project
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

    epi_cell_types <- c("CD8",
            "Mono_Macrophages",
            "CyclingImmune",
            "Plasma",
            "NK",
            "B Cells",
            "CD4",
            "Mast", "DC",
            "ILC","T Cells")

    # Can start with all cells filtered project and then subset
    # Define multiome subset to explore
    idxSample <- BiocGenerics::which((proj$Sample %in% unique(getCellColData(proj)$Sample)[order(unique(getCellColData(proj)$Sample))][27:71]) & 
        (proj$CellTypeRNA %in% epi_cell_types) &
        (proj$Sample %ni% c("B009-A-001" ,"B009-A-301")))
    cellsSample <- proj$cellNames[idxSample]

    proj_immune <- subsetArchRProject(
      ArchRProj = proj,
      cells = cellsSample,
      outputDirectory = new_project_save_name, dropCells = FALSE
    )

    saveArchRProject(ArchRProj = proj_immune, load = FALSE, overwrite = FALSE)
} else {
    proj_immune <- loadArchRProject(path = new_project_save_name)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 1) LSI Projection and Clustering
if (1 %in% execute_steps){
    proj_immune <- addIterativeLSI(
        ArchRProj = proj_immune,
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
    proj_immune <- addUMAP(
        ArchRProj = proj_immune, 
        reducedDims = paste("IterativeLSI", subscript, sep = ""), 
        name = paste("UMAP", subscript, sep = ""), 
        nNeighbors = 30, 
        minDist = 0.5, 
        metric = "cosine", force=TRUE
    )
    p1 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = "Sample", embedding = paste("UMAP", subscript, sep = ""))
    #p2 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = paste("Clusters", subscript, sep = ""), embedding = paste("UMAP", subscript, sep = ""))
    p3 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = "Donor", embedding = paste("UMAP", subscript, sep = ""))
    p4 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = "Location", embedding = paste("UMAP", subscript, sep = ""))
    plotPDF(p1,p3,p4, name = paste(paste("Plot-UMAP-Sample-Clusters-Donor-Location", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_immune, addDOC = FALSE, width = 5, height = 5)
    saveArchRProject(ArchRProj = proj_immune, load = FALSE, overwrite = FALSE)

    proj_immune <- addHarmony(
        ArchRProj = proj_immune,
        reducedDims = paste("IterativeLSI", subscript, sep = ""),
        name = paste("Harmony", subscript, sep = ""),
        groupBy = "Sample", force = TRUE
    )
    proj_immune <- addUMAP(
        ArchRProj = proj_immune, 
        reducedDims = paste("Harmony", subscript, sep = ""), 
        name = paste("UMAPHarmony", subscript, sep = ""), 
        nNeighbors = 30, 
        minDist = 0.5, 
        metric = "cosine", force = TRUE
    )
    proj_immune <- addClusters(
        input = proj_immune,
        reducedDims = paste("Harmony", subscript, sep = ""),
        method = "Seurat",
        name = paste("ClustersHarmony", subscript, sep = ""),
        resolution = 1.2, force=TRUE, nOutlier = 50, seed = 1, maxClusters = 30
    )
    p3 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = "Sample", embedding = paste("UMAPHarmony", subscript, sep = ""))
    p5 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = "Location", embedding = paste("UMAPHarmony", subscript, sep = ""))
    p6 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = "Donor", embedding = paste("UMAPHarmony", subscript, sep = ""))
    p7 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = paste("ClustersHarmony", subscript, sep = ""), embedding = paste("UMAPHarmony", subscript, sep = ""))
    plotPDF(p3,p5,p6,p7, name = paste(paste("Plot-UMAP2Harmony-Sample-Location-Donor-Clusters", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_immune, addDOC = FALSE, width = 5, height = 5)
    p6 <- plotEmbedding(ArchRProj = proj_immune, colorBy = "cellColData", name = "CellTypeRNA", embedding = paste("UMAPHarmony", subscript, sep = ""))
    plotPDF(p6, name = paste(paste("Plot-UMAP2Harmony-CellType", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj_immune, addDOC = FALSE, width = 5, height = 5)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 2) Add RNA data to multiome project
if (2 %in% execute_steps){
    colon <- readRDS(paste0(parent_directory, "scRNA/immune/decontx_norm_scale_harmony/clustered_full_colon_immune_proj_seurat_filtering_complete_counts_decontx_soup.rds"))

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

    rnaRowRanges <- readRDS("10x_RNA_row_ranges.rds")

    colon.data.full <- colon.data.full[rownames(colon.data.full) %in% names(rnaRowRanges),]
    library(stringr)
    colnames(colon.data.full) <- str_replace(colnames(colon.data.full), "_", "#")

    rowRangesSE <- rnaRowRanges[names(rnaRowRanges) %in% rownames(colon.data.full),]
    colon.data.full <- colon.data.full[names(rowRangesSE),] # make sure they have the same order

    seRNA<-SummarizedExperiment(assays=list(counts=colon.data.full), rowRanges= rowRangesSE) # create summarized experiment for adding to arrow

    proj_immune <- addGeneExpressionMatrix(input = proj_immune, seRNA = seRNA, force = TRUE)

    # save project
    saveArchRProject(ArchRProj = proj_immune, load = FALSE, overwrite = FALSE)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 3) Call Peaks
if (3 %in% execute_steps){
    proj_immune <- addGroupCoverages(ArchRProj = proj_immune, groupBy = "CellTypeRNA", force = TRUE)
    pathToMacs2 <- findMacs2()
    #Call Reproducible Peaks w/ Macs2
    proj_immune <- addReproduciblePeakSet(
        ArchRProj = proj_immune, groupBy = "CellTypeRNA", force = TRUE, 
        pathToMacs2 = pathToMacs2
    )
    #Add Peak Matrix
    proj_immune <- addPeakMatrix(ArchRProj = proj_immune, force = TRUE)
    saveArchRProject(ArchRProj = proj_immune, load = FALSE, overwrite = FALSE)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 4) Add deviations matrix
if (4 %in% execute_steps){
    proj_immune <- addBgdPeaks(proj_immune)
    motifPWMs <- readRDS("Vierstra-Human-Motifs.rds")
    proj_immune <- addMotifAnnotations(proj_immune, motifPWMs = motifPWMs, name = "Vierstra")
    proj_immune <- addDeviationsMatrix(proj_immune, peakAnnotation = "Vierstra", force = TRUE)
    saveArchRProject(ArchRProj = proj_immune, load = FALSE, overwrite = FALSE)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 5) Peak to gene linkages
if (5 %in% execute_steps){
    proj_immune <- addPeak2GeneLinks(
        ArchRProj = proj_immune,
        reducedDims = 
    )
}

