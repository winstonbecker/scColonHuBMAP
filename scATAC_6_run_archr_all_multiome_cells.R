# Create ArchR project for all multiome cells, and define a uniform peak set for all cell types

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
new_project_save_name = "HuBMAP_all_multiome"
subscript = "all_cells"
set.seed(1)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 1) Load and subset previously defined archr project

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

    idxSample <- BiocGenerics::which((proj$Sample %in% unique(getCellColData(proj)$Sample)[order(unique(getCellColData(proj)$Sample))][27:71]) & 
        (proj$Sample %ni% c("B009-A-001", "B009-A-301"))) # Note these two samples exhibited low diversity in the atac data and were excluded

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
# 1) Call Peaks
if (1 %in% execute_steps){
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

    cellLocationsEpi <- cellLocations[cellTypes$CellTypeRNA %in% epi_cell_types,, drop = FALSE]

    cellTypesEpi <- cellTypes[cellTypes$CellTypeRNA %in% epi_cell_types,, drop = FALSE]
    cellTypesImmStr <- cellTypes[cellTypes$CellTypeRNA %ni% epi_cell_types,, drop = FALSE]

    cellTypesEpi$CellTypeRNA <- paste0(cellLocationsEpi$Location, "-", cellTypesEpi$CellTypeRNA)
    cellTypesEpi[cellTypesEpi$CellTypeRNA == "Duodenum-Tuft",] <- "Proximal-jejunum-Tuft"
    cellTypesEpi[cellTypesEpi$CellTypeRNA == "Ileum-Tuft",] <- "Proximal-jejunum-Tuft"
    cellTypes <- rbind(cellTypesEpi, cellTypesImmStr)
    cellTypes[cellTypes$CellTypeRNA == "T Cells",] <- "CD4"

    proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste0(cellTypes$CellTypeRNA), cells = rownames(cellTypes), name = "CellTypeLocationPeakCalling", force = TRUE)

    set.seed(1)
    library(parallel)
    proj_stromal <- addGroupCoverages(ArchRProj = proj_stromal, groupBy = "CellTypeLocationPeakCalling", force = TRUE) #paste("Clusters", subscript, sep = "")
    pathToMacs2 <- findMacs2()

    #Call Reproducible Peaks w/ Macs2
    proj_stromal <- addReproduciblePeakSet(
        ArchRProj = proj_stromal, groupBy = "CellTypeLocationPeakCalling", force = TRUE, 
        pathToMacs2 = pathToMacs2, maxPeaks = 250000
    )

    #Add Peak Matrix
    proj_stromal <- addPeakMatrix(ArchRProj = proj_stromal, force = TRUE)

    saveArchRProject(ArchRProj = proj_stromal, load = FALSE, overwrite = FALSE)
}



