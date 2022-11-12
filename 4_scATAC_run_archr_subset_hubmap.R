# Script to create a normal intestine scATAC archr project

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

parent_directory <- "/hubmap_single_cell/"

# Set things up
.libPaths("./libraries/R_LIBS_4p1p2/")

# Load packages
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
`%notin%` <- Negate(`%in%`)

#Set/Create Working Directory to Folder
setwd(paste0(parent_directory, "scATAC/projects/"))

#Load Genome Annotations
addArchRGenome("hg38")

#Set Threads to be used
addArchRThreads()

# Things to set for subseting projects. Also find and replace proj name for subsample with desired project name (or similar)
subscript = "hubmap_all"

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
    "BATF","TNFRSF4", "FOXP3","CTLA4","LAIR2", # Tregs
    "FOLR2","FABP3","PLA2G2D", #Macrophages
    "FAP", "CBLN2", "SPOCK1", "ACSS3", # Fibroblast
    "SYT10", "SOSTDC1", "DES", "TAGLN", #Myofibroblasts
    "SELP", "ZNF385D", "FAM155A", "GALNT15", "MADCAM1", "CORT", #Post capillary venules
    "COX4I2", "KCNJ8", "HIGD1B", "RGS5", "NOTCH3", "HEYL", "FAM162B", #Pericytes
    "FAM110D", "INHBB", "NPR1", "NOVA2", "GPIHBP1", "SOX17", #endothelial
    "S100A1", # nerves
    "RSPO3", "CCL11", "WNT5B", "BMP4", "CHI3L1", "ACTA2", "WNT2B", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", 
    "MADCAM1", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A",
    "DCLK1", #Tuft
    "HTR3C", "HTR3E", "B4GALNT4", # Tuft
    "KLK1","ITLN1","WFDC2","CLCA1", # Immature Goblet
    "MUC2", "TFF1", "FCGBP","TBX10", # Goblet
    "FAP", # Fibroblast
    "CA1", # E.Immature_Enterocytes
    "RAB6B", #Enterocytes
    "CRYBA2","SCGN", #Enteroendocrine
    "CA2", "SI", # absorptive
    "SOX9", "CD34", #progenitor
    "MUC1", "KRT1", # General epithelial
    "LYZ", "DEFA5", # Paneth
    "GP2", "KRT7", # M cells
    "NTRK2","CCL23", # M cells
    "BEST4", "CA7","OTOP2","OTOP3", # Best4+ enterocytes
    "S100A1", # nerves
    "SMOC2", "RGMB", "LGR5", "ASCL2" #stem
)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 1) Load Project
proj <- loadArchRProject("all_hubmap_cells")
proj <- filterDoublets(proj, filterRatio = 1.2)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 2) Add metadata
metadata <- read.table(paste0(parent_directory, "hubmap_metadata_atac_and_rna_final.csv"), header = TRUE, sep = ",", stringsAsFactors=FALSE)
for (j in 2:dim(metadata)[2]){
  # initialize list
  cellsNamesToAdd <- c()
  annotationToAdd <- c()
  for (i in 1:dim(metadata)[1]){
    idxSample <- BiocGenerics::which(getCellColData(proj, "Sample") %in% metadata[i,"Sample"])
    cellsSample <- proj$cellNames[idxSample[["Sample"]]]
    cellsNamesToAdd <- append(cellsNamesToAdd, cellsSample)
    annotationToAdd <- append(annotationToAdd, rep(metadata[i,j], length(cellsSample)))
  }
  proj <- addCellColData(ArchRProj = proj, data = paste0(annotationToAdd), cells = paste0(cellsNamesToAdd), name = colnames(metadata)[j], force = TRUE)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 2) load scrna data and filter out multiome cells that were filtered out in rna analysis
immune <- readRDS(paste0(parent_directory, "scRNA/immune/decontx_norm_scale_harmony/diet_clustered_full_colon_immune_proj_seurat_filtering_complete.rds"))
stromal <- readRDS(paste0(parent_directory, "scRNA/stromal/decontx_norm_scale_harmony/diet_stromal_all_samples_clustered_filtered.rds"))
colon <- readRDS(paste0(parent_directory, "scRNA/epithelial/colon/decontx_sctransform_cca/colon_clustered.rds"))
ileum <- readRDS(paste0(parent_directory, "scRNA/epithelial/ileum/decontx_sctransform_cca/ileum.rds"))
jejunum <- readRDS(paste0(parent_directory, "scRNA/epithelial/jejunum/decontx_sctransform_cca/jejunum_clustered.rds"))
duodenum <- readRDS(paste0(parent_directory, "scRNA/epithelial/duodenum/decontx_sctransform_cca/duodenum_clustered.rds"))
enteroendocrine <- readRDS(paste0(parent_directory, "scRNA/epithelial/enteroendocrine/decontx_sctransform_cca/clustered_annotated_enteroendocrine_cells.rds"))
secretory_special <- readRDS(paste0(parent_directory, "scRNA/epithelial/specialized_secretory/decontx_sctransform_cca/specialized_secretory_clustered.rds"))

# get the cell types
celltypes <- rbind(immune@meta.data[,"CellType", drop = FALSE],
    stromal@meta.data[,"CellType", drop = FALSE],
    ileum@meta.data[,"CellType", drop = FALSE],
    colon@meta.data[,"CellType", drop = FALSE],
    jejunum@meta.data[,"CellType", drop = FALSE],
    duodenum@meta.data[,"CellType", drop = FALSE])

# keep the more specific cell types
celltypes_sp <- rbind(enteroendocrine@meta.data[,"CellType", drop = FALSE],
  secretory_special@meta.data[,"CellType", drop = FALSE])
celltypes_full <- rbind(celltypes[!(rownames(celltypes) %in% rownames(celltypes_sp)),, drop = FALSE],
  celltypes_sp)
write.table(celltypes_full, "scrna_cell_types.tsv")

# get the samples, multiome samples, and non multiome samples
samples <- unique(c(paste0(immune$orig.ident),
    paste0(stromal$orig.ident),
    paste0(colon$orig.ident),
    paste0(ileum$orig.ident),
    paste0(jejunum$orig.ident),
    paste0(duodenum$orig.ident)))

non_multiome <- c(samples[grepl("B001", samples)],
  samples[grepl("B004", samples)],
  samples[grepl("B005", samples)], "B001-A-302", "B004-A-004-R2")

multiome <- c(samples[grepl("B006", samples)],
  samples[grepl("B008", samples)],
  samples[grepl("B009", samples)],
  samples[grepl("B010", samples)],
  samples[grepl("B011", samples)],
  samples[grepl("B012", samples)])

# create a list of all the cells
all_cells <- c(colnames(immune), colnames(stromal), colnames(colon), colnames(ileum), colnames(jejunum), colnames(duodenum))
library(stringr)
all_cells <- str_replace(all_cells, "_", "#")

# create a list of just the multiome cells
all_cells <- all_cells[(grepl(paste(multiome, collapse="|"), all_cells))]

# keep all the nonmultiome cells and the multiome cells with rna annotations
keep <- c(getCellNames(proj)[(grepl(paste(non_multiome, collapse="|"), getCellNames(proj)))],
  getCellNames(proj)[getCellNames(proj) %in% all_cells])
table(substring(keep, 1, nchar(keep)-19))

# save a new filtered project
new_project_save_name <- "all_hubmap_cells_filtered"
proj_filtered <- subsetArchRProject(
  ArchRProj = proj,
  cells = keep,
  outputDirectory = new_project_save_name, dropCells = FALSE
)
proj <- loadArchRProject("all_hubmap_cells_filtered")


################################################################################################################
#..............................................................................................................#
################################################################################################################
# 3) QC Violin Plots
p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.8,
    addBoxPlot = TRUE
)
p3 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "nFrags",
    plotAs = "violin",
    alpha = 0.8,
    addBoxPlot = TRUE
)
plotPDF(p2,p3, name = "QC-Sample-Statistics-TSS-Sample-Violin-Doublets-Filtered.pdf", ArchRProj = proj, addDOC = FALSE, width = 20, height = 15)

pal <- paletteDiscrete(unique(getCellColData(proj)$Sample))
palMap <- c("Ascending" = "#208A42", "Descending" = "#272E6A", "Duodenum" = "#D51F26", "Ileum" = "#FEE500", 
  "Mid-jejunum" = "#F8AD16", "Proximal-jejunum" = "#F47D2B", "Sigmoid" = "#89288F", "Transverse" = "#438496")
for (i in names(pal)){
  pal[i] <- paste0(palMap[metadata[metadata$Sample == i,]$Location])
}
p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    name = "TSSEnrichment",
    colorBy = "cellColData", 
    plotAs = "violin",
    alpha = 1,
    addBoxPlot = FALSE, pal = pal)
#p2 <- p2+geom_boxplot(outlier.shape = NA, alpha = 1)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))

p2 <- p2+geom_boxplot(outlier.shape = NA, alpha = 1)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
      scale_x_discrete(labels= paste0(names(table(proj$Sample)), " (n=", table(proj$Sample), ")"))

plotPDF(p2, name = "Temp-QC-Sample-Statistics-TSS-Sample-Violin-Doublets-Filtered-with-N.pdf", ArchRProj = proj, addDOC = FALSE, width = 20, height = 15)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 4) Dimensionality reduciton and clustering
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = paste("IterativeLSI", subscript, sep = ""), 
    iterations = 3, 
    clusterParams = list(
        resolution = c(0.1, 0.2), 
        sampleCells = 20000, 
        n.start = 10
    ), 
    varFeatures = 15000, sampleCellsPre = 25000,
    dimsToUse = 1:25, force = TRUE
)

# # alternatively, read in already defined iterative LSI
# IterativeLSIhubmap_all <- readRDS("all_hubmap_cells_filtered__IterativeLSIhubmap_all.rds")
# Harmonyhubmap_all <- readRDS("all_hubmap_cells_filtered__Harmonyhubmap_all.rds")
# clusters <- readRDS("initial_clusters.rds")

# proj@reducedDims[["IterativeLSIhubmap_all"]] <- IterativeLSIhubmap_all
# proj@reducedDims[["Harmonyhubmap_all"]] <- Harmonyhubmap_all
# proj <- addCellColData(ArchRProj = proj, data = paste0(clusters$Clustershubmap_all), cells = paste0(rownames(clusters)), name = "Clustershubmap_all", force = TRUE)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = paste("IterativeLSI", subscript, sep = ""), 
    name = paste("UMAP", subscript, sep = ""), 
    nNeighbors = 30, 
    minDist = 0.4, 
    metric = "cosine", force=TRUE
)
proj <- addClusters(
    input = proj,
    reducedDims = paste("IterativeLSI", subscript, sep = ""),
    method = "Seurat",
    name = paste("Clusters", subscript, sep = ""),
    resolution = 1.5, force=TRUE, nOutlier = 20, seed = 1, sampleCells = 40000, maxClusters = 40
)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = paste("UMAP", subscript, sep = ""))
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Donor", embedding = paste("UMAP", subscript, sep = ""))

pal <- c("Ascending" = "#208A42", "Descending" = "#272E6A", "Duodenum" = "#D51F26", "Ileum" = "#FEE500", 
  "Mid-jejunum" = "#F8AD16", "Proximal-jejunum" = "#F47D2B", "Sigmoid" = "#89288F", "Transverse" = "#438496")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Location", embedding = paste("UMAP", subscript, sep = ""))#, pal = pal)
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = paste("Clusters", subscript, sep = ""), embedding = paste("UMAP", subscript, sep = ""))
plotPDF(p1,p2,p3,p4, name = paste(paste("Plot-UMAP-Sample-Donor-Location-Clusters", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DoubletScore", embedding = paste("UMAP", subscript, sep = ""))#, pal = pal)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = paste("UMAP", subscript, sep = ""))#, pal = pal)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "nFrags", embedding = paste("UMAP", subscript, sep = ""))#, pal = pal)
plotPDF(p1,p2,p3, name = paste(paste("Plot-UMAP-QC-Stats", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj, outputDirectory = "all_hubmap_cells_filtered", load = FALSE, overwrite = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 5) Plot Gene Scores

proj <- addImputeWeights(proj, reducedDims = paste0("IterativeLSI", subscript))
p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = paste("UMAP", subscript, sep = ""),
    imputeWeights = getImputeWeights(proj)
)
plotPDF(plotList = p, 
    name = paste(paste("Plot-UMAP-IterativeLSI-Marker-Genes-W-Imputation", subscript, sep = "-"), ".pdf", sep = ""), 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Harmony, clustering, and cell type labeling
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = paste("IterativeLSI", subscript, sep = ""),
    name = paste("Harmony", subscript, sep = ""),
    groupBy = "Donor", force = TRUE
)
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = paste("Harmony", subscript, sep = ""), 
    name = paste("UMAPHarmony", subscript, sep = ""), 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine", force = TRUE
)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = paste("UMAPHarmony", subscript, sep = ""))
p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Location", embedding = paste("UMAPHarmony", subscript, sep = ""))
p6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Donor", embedding = paste("UMAPHarmony", subscript, sep = ""))
plotPDF(p3,p5,p6, name = paste(paste("Plot-UMAP2Harmony-Sample-Location-Donor", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj, outputDirectory = "all_hubmap_cells_filtered", load = FALSE, overwrite = FALSE)


################################################################################################################
#..............................................................................................................#
################################################################################################################
# Plot cell types
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

p6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "CellTypeRNA", embedding = paste("UMAPHarmony", subscript, sep = ""))
plotPDF(p6, name = paste(paste("Plot-UMAP2Harmony-CellTypeRNA", subscript, sep = "-"), ".pdf", sep = ""), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

cell_types_scrna <- read.table("scrna_cell_types.tsv")
samples <- unique(substr(rownames(cell_types_scrna), 1, nchar(rownames(cell_types_scrna))-19))

# Save
saveArchRProject(ArchRProj = proj, outputDirectory = "all_hubmap_cells_filtered", load = FALSE, overwrite = FALSE)



