# Script to identify possible regulatory transcription factors in the single cell hubmap data
# WRB 2022

# Set things up
.libPaths("./libraries/R_LIBS_4p1p2/")

# Load packages
library(ArchR)
library(Seurat)
library(parallel)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
`%notin%` <- Negate(`%in%`)

#Load Hidden Helper Functions
fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
for (i in seq_along(fn)) {
    tryCatch({
        eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
    }, error = function(x) {
    })
}

#Set/Create Working Directory to Folder
parent_directory <- "/hubmap_single_cell/"
setwd(paste0(parent_directory, "scATAC/projects/"))

#Load Genome Annotations
addArchRGenome("hg38")

#Set Threads to be used
addArchRThreads()

# Set correlation cutoff for TF regulators
corr_cutoff = 0.5
corr_savename <- str_replace(paste0(corr_cutoff), "0.", "0p")

# Load previously defined archr projects (each project should have gene expression matrix and vierstra motif matrix)
# for epithelial
proj_ileum <- loadArchRProject(path = "HuBMAP_ileum_epithelial_cells_multiome")
proj_colon <- loadArchRProject(path = "HuBMAP_colon_epithelial_cells_multiome")
proj_duodenum <- loadArchRProject(path = "HuBMAP_duodenum_epithelial_cells_multiome")
proj_jejunum <- loadArchRProject(path = "HuBMAP_jejunum_epithelial_cells_multiome")
proj_epithelial <- loadArchRProject(path = "HuBMAP_epithelial_cells_multiome_final")
# for immune
proj_immune <- loadArchRProject(path = "HuBMAP_immune_cells_multiome")
# for stromal
proj_stromal <- loadArchRProject(path = "HuBMAP_stromal_cells_multiome")

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 2) Define Functions
get_TF_regulators <- function(proj_stromal, saveName, useMatrix = "VierstraMatrix", groupBy = "CellTypeRNA", corr = 0.5, padj_cut = 0.01, delta_quant = 0.7){
  seGroupMotif <- getGroupSE(ArchRProj = proj_stromal, useMatrix = useMatrix, groupBy = groupBy)
  seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
  rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
    rowMaxs(assay(seZ) - assay(seZ)[,x])
  }) %>% Reduce("cbind", .) %>% rowMaxs

  corGIM_MM <- correlateMatrices(
      ArchRProj = proj_stromal,
      useMatrix1 = "GeneExpressionMatrix",
      useMatrix2 = "VierstraMatrix",
      reducedDims = paste0("Harmony", saveName)
  )

  corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$VierstraMatrix_name, rowData(seZ)$name), "maxDelta"]

  corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
  corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,paste0(useMatrix, "_name")]))), ]
  corGIM_MM$TFRegulator <- "NO"
  corGIM_MM$TFRegulator[which(corGIM_MM$cor > corr & corGIM_MM$padj < padj_cut & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, delta_quant))] <- "YES"
  sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
  p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator, label = GeneExpressionMatrix_name)) +
    geom_point() +
    geom_text(data = data.frame(subset(corGIM_MM, corGIM_MM$cor > corr & corGIM_MM$padj < padj_cut & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, delta_quant))), hjust=-0.2, vjust=0.2)+
    theme_ArchR() +
    geom_vline(xintercept = 0, lty = "dashed") + 
    scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
    xlab("Correlation To Gene Expression") +
    ylab("Max TF Motif Delta")+
    xlim(-1.3, 1.3) +
    scale_y_continuous(
      expand = c(0,0), 
      limits = c(0, max(corGIM_MM$maxDelta)*1.1)
    )
  plotPDF(p, name = paste(paste(saveName, "positive_TF_regulators_gene_integration_matrix", sep = "-"), ".pdf", sep = ""), width = 8, height = 8, ArchRProj = proj_stromal, addDOC = FALSE)

  motifs <- unique(sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1]))

  p <- plotEmbedding(
      ArchRProj = proj_stromal, 
      colorBy = "GeneExpressionMatrix", 
      name = motifs,
      pal = paletteContinuous("blueYellow"),
      embedding = paste("UMAPHarmony", saveName, sep = ""),
      imputeWeights = getImputeWeights(proj_stromal)
  )
  plotPDF(p, name = paste0("Gene-Expression-TF-regulators-", saveName, ".pdf"), width = 5, height = 5,  ArchRProj = proj_stromal, addDOC = FALSE)

  return(corGIM_MM[corGIM_MM$TFRegulator=="YES",1:2])
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Epithelial Cells

# ID regulators for each section individually
motifs_duodenum <- get_TF_regulators(proj_stromal = proj_duodenum, saveName = "epithelial_duodenum", corr = corr_cutoff)
motifs_jejunum <- get_TF_regulators(proj_stromal = proj_jejunum, saveName = "epithelial_jejunum", corr = corr_cutoff)
motifs_ileum <- get_TF_regulators(proj_stromal = proj_ileum, saveName = "epithelial_ileum", corr = corr_cutoff)
motifs_colon <- get_TF_regulators(proj_stromal = proj_colon, saveName = "epithelial_colon", corr = corr_cutoff)
motifs_epithelial <- get_TF_regulators(proj_stromal = proj_epithelial, saveName = "hubmap_all", corr = corr_cutoff)

# Create combined list of regulators in any region (only keeping one motif for the same TF)
all_motifs_both_names <- rbind(motifs_duodenum, motifs_jejunum, motifs_ileum, motifs_colon, motifs_epithelial)
all_motifs_both_names <- all_motifs_both_names[! duplicated(all_motifs_both_names$GeneExpressionMatrix_name),]
all_motifs_both_names <- all_motifs_both_names[all_motifs_both_names$GeneExpressionMatrix_name %ni% c("HNF1A-AS1", "RFX3-AS1", "MEIS1-AS2", "HNF4A-AS1", "MSC-AS1", "SOX2-OT", "TGIF2-RAB5IF", "SOX9-AS1", "HNF4A-AS1", "ELF3-AS1", "GATA6-AS1"),]

# Get gene expression matrix for each region
GIM_duodenum <- getGroupSE(ArchRProj = proj_duodenum, useMatrix = "GeneExpressionMatrix", groupBy = "CellTypeRNA")
GIM_jejunum <- getGroupSE(ArchRProj = proj_jejunum, useMatrix = "GeneExpressionMatrix", groupBy = "CellTypeRNA")
GIM_ileum <- getGroupSE(ArchRProj = proj_ileum, useMatrix = "GeneExpressionMatrix", groupBy = "CellTypeRNA")
GIM_colon <- getGroupSE(ArchRProj = proj_colon, useMatrix = "GeneExpressionMatrix", groupBy = "CellTypeRNA")

# Get motif matrix for everything together - note this is important because you want the version where the deviations were all calculated together
# When doing this, we combine serveral of the specific cell types that have relatively few cells represented across the colon
locations <- getCellColData(proj_epithelial)$Location
locations[locations == "Ascending"] <- "colon"
locations[locations == "Transverse"] <- "colon"
locations[locations == "Descending"] <- "colon"
locations[locations == "Sigmoid"] <- "colon"
locations[locations == "Duodenum"] <- "duodenum"
locations[locations == "Proximal-jejunum"] <- "jejunum"
locations[locations == "Mid-jejunum"] <- "jejunum"
locations[locations == "Ileum"] <- "ileum"
cellTypeLocation <- paste(getCellColData(proj_epithelial)$CellTypeRNA, locations, sep = "_")
cellTypeLocation[cellTypeLocation == "D Cells_jejunum"] <- "D Cells"
cellTypeLocation[cellTypeLocation == "D Cells_ileum"] <- "D Cells"
cellTypeLocation[cellTypeLocation == "D Cells_colon"] <- "D Cells"
cellTypeLocation[cellTypeLocation == "D Cells_duodenum"] <- "D Cells"
cellTypeLocation[cellTypeLocation == "S Cells_jejunum"] <- "S Cells"
cellTypeLocation[cellTypeLocation == "S Cells_ileum"] <- "S Cells"
cellTypeLocation[cellTypeLocation == "S Cells_colon"] <- "S Cells"
cellTypeLocation[cellTypeLocation == "S Cells_duodenum"] <- "S Cells"
cellTypeLocation[cellTypeLocation == "I Cells_jejunum"] <- "I Cells"
cellTypeLocation[cellTypeLocation == "I Cells_ileum"] <- "I Cells"
cellTypeLocation[cellTypeLocation == "I Cells_colon"] <- "I Cells"
cellTypeLocation[cellTypeLocation == "I Cells_duodenum"] <- "I Cells"
cellTypeLocation[cellTypeLocation == "K Cells_jejunum"] <- "K Cells"
cellTypeLocation[cellTypeLocation == "K Cells_ileum"] <- "K Cells"
cellTypeLocation[cellTypeLocation == "K Cells_colon"] <- "K Cells"
cellTypeLocation[cellTypeLocation == "K Cells_duodenum"] <- "K Cells"
cellTypeLocation[cellTypeLocation == "L Cells_jejunum"] <- "L Cells"
cellTypeLocation[cellTypeLocation == "L Cells_ileum"] <- "L Cells"
cellTypeLocation[cellTypeLocation == "L Cells_colon"] <- "L Cells"
cellTypeLocation[cellTypeLocation == "L Cells_duodenum"] <- "L Cells"
cellTypeLocation[cellTypeLocation == "Mo Cells_jejunum"] <- "Mo Cells"
cellTypeLocation[cellTypeLocation == "Mo Cells_ileum"] <- "Mo Cells"
cellTypeLocation[cellTypeLocation == "Mo Cells_colon"] <- "Mo Cells"
cellTypeLocation[cellTypeLocation == "Mo Cells_duodenum"] <- "Mo Cells"
cellTypeLocation[cellTypeLocation == "NEUROG3high_jejunum"] <- "NEUROG3high"
cellTypeLocation[cellTypeLocation == "NEUROG3high_ileum"] <- "NEUROG3high"
cellTypeLocation[cellTypeLocation == "NEUROG3high_colon"] <- "NEUROG3high"
cellTypeLocation[cellTypeLocation == "NEUROG3high_duodenum"] <- "NEUROG3high"
cellTypeLocation[cellTypeLocation == "EnteroendocrineUn 1_jejunum"] <- "EnteroendocrineUn 1"
cellTypeLocation[cellTypeLocation == "EnteroendocrineUn 1_ileum"] <- "EnteroendocrineUn 1"
cellTypeLocation[cellTypeLocation == "EnteroendocrineUn 1_colon"] <- "EnteroendocrineUn 1"
cellTypeLocation[cellTypeLocation == "EnteroendocrineUn 1_duodenum"] <- "EnteroendocrineUn 1"
cellTypeLocation[cellTypeLocation == "EnteroendocrineUn_jejunum"] <- "EnteroendocrineUn"
cellTypeLocation[cellTypeLocation == "EnteroendocrineUn_ileum"] <- "EnteroendocrineUn"
cellTypeLocation[cellTypeLocation == "EnteroendocrineUn_colon"] <- "EnteroendocrineUn"
cellTypeLocation[cellTypeLocation == "EnteroendocrineUn_duodenum"] <- "EnteroendocrineUn"
cellTypeLocation[cellTypeLocation == "Enterochromaffin_jejunum"] <- "Enterochromaffin"
cellTypeLocation[cellTypeLocation == "Enterochromaffin_ileum"] <- "Enterochromaffin"
cellTypeLocation[cellTypeLocation == "Enterochromaffin_colon"] <- "Enterochromaffin"
cellTypeLocation[cellTypeLocation == "Enterochromaffin_duodenum"] <- "Enterochromaffin"

#proj_epithelial <- addCellColData(proj_epithelial, data = paste(getCellColData(proj_epithelial)$CellTypeRNA, locations, sep = "_"), name = "cellTypeLocationSimple", cells = rownames(getCellColData(proj_epithelial)))
proj_epithelial <- addCellColData(proj_epithelial, data = cellTypeLocation, name = "cellTypeLocationSimple", cells = rownames(getCellColData(proj_epithelial)), force = TRUE)
MM_epithelail <- getGroupSE(ArchRProj = proj_epithelial, useMatrix = "VierstraMatrix", groupBy = "cellTypeLocationSimple")
GIM_epithelail <- getGroupSE(ArchRProj = proj_epithelial, useMatrix = "GeneExpressionMatrix", groupBy = "cellTypeLocationSimple")

# Subset for things that are regulators
GIM_duodenum_sub <- GIM_duodenum[rowData(GIM_duodenum)$name %in% all_motifs_both_names$GeneExpressionMatrix_name,]
GIM_jejunum_sub <- GIM_jejunum[rowData(GIM_jejunum)$name %in% all_motifs_both_names$GeneExpressionMatrix_name,]
GIM_ileum_sub <- GIM_ileum[rowData(GIM_ileum)$name %in% all_motifs_both_names$GeneExpressionMatrix_name,]
GIM_colon_sub <- GIM_colon[rowData(GIM_colon)$name %in% all_motifs_both_names$GeneExpressionMatrix_name,]
GIM_epithelail_sub <- GIM_epithelail[rowData(GIM_epithelail)$name %in% all_motifs_both_names$GeneExpressionMatrix_name,]
MM_epithelail_sub <- MM_epithelail[rowData(MM_epithelail)$name %in% all_motifs_both_names$VierstraMatrix_name,]

# Make a combined expression matrix with all cell types for plotting
# test <- cbind(assays(GIM_duodenum_sub)$GeneExpressionMatrix, assays(GIM_jejunum_sub)$GeneExpressionMatrix, assays(GIM_ileum_sub)$GeneExpressionMatrix, assays(GIM_colon_sub)$GeneExpressionMatrix)
# rownames(test) <- rowData(GIM_duodenum_sub)$name
# colnames(test) <- c(paste0(colnames(GIM_duodenum_sub), "_duodenum"),paste0(colnames(GIM_jejunum_sub), "_jejunum"),
#     paste0(colnames(GIM_ileum_sub), "_ileum"),paste0(colnames(GIM_colon_sub), "_colon"))
# test <- test[all_motifs_both_names$GeneExpressionMatrix_name,]
test <- assays(GIM_epithelail_sub)$GeneExpressionMatrix
rownames(test) <- rowData(GIM_epithelail_sub)$name
test <- test[all_motifs_both_names$GeneExpressionMatrix_name,]

# Make a combined motif matrix with all cell types for plotting (just want the z-scores)
test_MM <- assays(MM_epithelail_sub)$VierstraMatrix[(length(rowData(MM_epithelail_sub)$name)/2+1):length(rowData(MM_epithelail_sub)$name),]
rownames(test_MM) <- rowData(MM_epithelail_sub)$name[(length(rowData(MM_epithelail_sub)$name)/2+1):length(rowData(MM_epithelail_sub)$name)]
test_MM <- test_MM[all_motifs_both_names$VierstraMatrix_name,]

# identify the order based on clustering both the motif matrix and gene score matrix at the same time
temp <- cbind(.rowZscores(test), .rowZscores(test_MM))
order <- hclust(dist(as.matrix(temp)))$order

# plot heatmaps
new_colnames <- c("Goblet_duodenum", "Goblet_jejunum", "Goblet_ileum", "Goblet_colon",
      "Immature Goblet_duodenum", "Immature Goblet_jejunum", "Immature Goblet_ileum", "Immature Goblet_colon",
      "Enterochromaffin",
      "EnteroendocrineUn",
      "EnteroendocrineUn 1",
      "D Cells",
      "I Cells",
      "K Cells",
      "L Cells",
      "Mo Cells",
      "NEUROG3high",
      "S Cells",
      "Exocrine_duodenum", "Secretory Specialized MUC5B+_duodenum", "Secretory Specialized MUC6+_duodenum", "Unknown_duodenum",
      "Paneth_duodenum", "Paneth_jejunum", "Paneth_ileum", 
      "Tuft_duodenum", "Tuft_jejunum", "Tuft_ileum","Tuft_colon",
      "Stem_duodenum", "Stem_jejunum", "Stem_ileum", "Stem_colon", 
      "CyclingTA 1_duodenum", "CyclingTA 1_jejunum", "CyclingTA 1_ileum", "CyclingTA 1_colon", "CyclingTA 2_colon", 
      "TA2_duodenum", "TA2_jejunum", "TA2_ileum", "TA2_colon", 
      "TA1_duodenum", "TA1_jejunum", "TA1_ileum", "TA1_colon", 
      "Immature Enterocytes_duodenum", "Immature Enterocytes_jejunum", "Immature Enterocytes_ileum", "Immature Enterocytes_colon", 
      "Enterocytes_duodenum", "Enterocytes_jejunum", "Enterocytes_ileum", "Enterocytes_colon",
      "Best4+ Enterocytes_duodenum", "Best4+ Enterocytes_jejunum", "Best4+ Enterocytes_ileum", "Best4+ Enterocytes_colon", 
      "Epithelial_duodenum", "Epithelial_jejunum", "Epithelial_ileum")#, "Epithelial_colon")

# subset to only keep things with at least 25 cells per group
new_colnames <- new_colnames[new_colnames %in% names(table(proj_epithelial$cellTypeLocationSimple)[table(proj_epithelial$cellTypeLocationSimple)>25])]

test <- test[,new_colnames]
test <- .rowZscores(test)
ht1 <- .ArchRHeatmap(
    mat = test[order,],
    limits = c(-2,2),
    color = paletteContinuous(set = "blueYellow", n = 200), 
    clusterCols = FALSE, 
    clusterRows = FALSE,
    labelRows = TRUE,
    labelCols = TRUE,
    customColLabel = NULL,
    showRowDendrogram = FALSE,
    showColDendrogram = FALSE,
    draw = FALSE,
    #split = clustering$cluster, #scale = TRUE,
    name = paste0(nrow(test), " features\n")
)
plotPDF(ht1, name = paste0("Gene-Integration-TF-regulators-heatmap-", "all_sections", corr_savename, "corrcut-ee-subs-combined", ".pdf"), width = 14, height = 15,  ArchRProj = proj_epithelial, addDOC = FALSE)

test_MM <- test_MM[,new_colnames]
test_MM <- .rowZscores(test_MM)
ht1 <- .ArchRHeatmap(
    mat = test_MM[order,],
    limits = c(-2,2),
    color = paletteContinuous(set = "solarExtra", n = 200), 
    clusterCols = FALSE, 
    clusterRows = FALSE,
    labelRows = TRUE,
    labelCols = TRUE,
    customColLabel = NULL,
    showRowDendrogram = FALSE,
    showColDendrogram = FALSE,
    draw = FALSE,
    #split = clustering$cluster, #scale = TRUE,
    name = paste0(nrow(test_MM), " features\n")
)
plotPDF(ht1, name = paste0("Motif-Matrix-TF-regulators-heatmap-", "all_sections", corr_savename, "corrcut-ee-subs-combined", ".pdf"), width = 14, height = 15,  ArchRProj = proj_epithelial, addDOC = FALSE)

write.csv(test[order,], paste0(parent_directory, "scATAC/projects/HuBMAP_epithelial_cells_multiome_final/Plots/Gene-Expression-Matrix-TF-regulators-heatmap-all_sections0p5corrcut-ee-subs-combined-r.csv"))
write.csv(test_MM[order,], paste0(parent_directory, "scATAC/projects/HuBMAP_epithelial_cells_multiome_final/Plots/Motif-Matrix-TF-regulators-heatmap-all_sections0p5corrcut-ee-subs-combined-r.csv"))


# Plot Footprints Example
motifPositions <- getPositions(proj_epithelial)
motifs <- c("RUNX2")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
seFoot <- getFootprints(
  ArchRProj = proj_epithelial, 
  positions = motifPositions[markerMotifs], 
  groupBy = "CellTypeLocation"
)
seFootsub <- seFoot[,grepl("Proximal-jejunum-Tuft", colnames(seFoot)) | grepl("Proximal-jejunum-Enterocytes", colnames(seFoot))]
p <- plotFootprints(
  seFoot = seFootsub,
  ArchRProj = proj_epithelial, 
  normMethod = "Divide",
  plotName = "Footprints",
  addDOC = FALSE,
  smoothWindow = 5, plot = FALSE
)
plotPDF(p, name = paste0(motifs, "_motif_footprinting_all_stem", ".pdf"), width = 14, height = 15,  ArchRProj = proj_epithelial, addDOC = FALSE)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Immune Cells
cellTypes <- getCellColData(proj_immune)[,c("CellTypeRNA"), drop = FALSE]
cellLocations <- getCellColData(proj_immune)[,c("Location"), drop = FALSE]
cellLocations[cellLocations$Location == "Transverse",] <- "Colon"
cellLocations[cellLocations$Location == "Sigmoid",] <- "Colon"
cellLocations[cellLocations$Location == "Descending",] <- "Colon"
cellLocations[cellLocations$Location == "Ascending",] <- "Colon"
cellLocations[cellLocations$Location == "Mid-jejunum",] <- "Proximal-jejunum"

cellTypes$CellTypeRNA <- paste0(cellLocations$Location, "-", cellTypes$CellTypeRNA)
cellTypes[cellTypes$CellTypeRNA == "Duodenum-T Cells",] <- "T Cells"
cellTypes[cellTypes$CellTypeRNA == "Proximal-jejunum-T Cells",] <- "T Cells"
cellTypes[cellTypes$CellTypeRNA == "Ileum-T Cells",] <- "T Cells"
cellTypes[cellTypes$CellTypeRNA == "Colon-T Cells",] <- "T Cells"

cellTypes[cellTypes$CellTypeRNA == "Duodenum-ILC",] <- "ILC"
cellTypes[cellTypes$CellTypeRNA == "Proximal-jejunum-ILC",] <- "ILC"
cellTypes[cellTypes$CellTypeRNA == "Ileum-ILC",] <- "ILC"
cellTypes[cellTypes$CellTypeRNA == "Colon-ILC",] <- "ILC"

cellTypes[cellTypes$CellTypeRNA == "Duodenum-DC",] <- "DC"
cellTypes[cellTypes$CellTypeRNA == "Proximal-jejunum-DC",] <- "DC"
cellTypes[cellTypes$CellTypeRNA == "Ileum-DC",] <- "DC"
cellTypes[cellTypes$CellTypeRNA == "Colon-DC",] <- "DC"

proj_immune <- addCellColData(ArchRProj = proj_immune, data = paste0(cellTypes$CellTypeRNA), cells = rownames(cellTypes), name = "CellTypeLocation", force = TRUE)

motifs_immune <- get_TF_regulators(proj_stromal = proj_immune, saveName = "immune", groupBy = "CellTypeLocation", corr = corr_cutoff)
motifs_immune <- motifs_immune[! duplicated(motifs_immune$GeneExpressionMatrix_name),]
# dont plot antisense rnas correlated to motif activity
motifs_immune <- motifs_immune[motifs_immune$GeneExpressionMatrix_name %ni% c("HNF1A-AS1", "RFX3-AS1", "MEIS1-AS2", "HNF4A-AS1", "MSC-AS1", "SOX2-OT", "TGIF2-RAB5IF","ETV5-AS1", "MEF2C-AS1", "RORA-AS1", "CEBPB-AS1"),]

MM_epithelail <- getGroupSE(ArchRProj = proj_immune, useMatrix = "VierstraMatrix", groupBy = "CellTypeLocation")
GIM_epithelail <- getGroupSE(ArchRProj = proj_immune, useMatrix = "GeneExpressionMatrix", groupBy = "CellTypeLocation")

MM_epithelail_sub <- MM_epithelail[rowData(MM_epithelail)$name %in% motifs_immune$VierstraMatrix_name,]
GIM_epithelail_sub <- GIM_epithelail[rowData(GIM_epithelail)$name %in% motifs_immune$GeneExpressionMatrix_name,]


# Make a combined matrix with all cell types for plotting
test <- assays(GIM_epithelail_sub)$GeneExpressionMatrix
rownames(test) <- rowData(GIM_epithelail_sub)$name
colnames(test) <- colnames(GIM_epithelail_sub)
test <- test[motifs_immune$GeneExpressionMatrix_name,]

# get deviation z scores, which is the second half of the matrix
test_MM <- assays(MM_epithelail_sub)$VierstraMatrix[(length(rowData(MM_epithelail_sub)$name)/2+1):length(rowData(MM_epithelail_sub)$name),]
rownames(test_MM) <- rowData(MM_epithelail_sub)$name[(length(rowData(MM_epithelail_sub)$name)/2+1):length(rowData(MM_epithelail_sub)$name)]
test_MM <- test_MM[motifs_immune$VierstraMatrix_name,]

# identify the order based on both the motif matrix and gene score matrix at the same time
temp <- cbind(.rowZscores(test), .rowZscores(test_MM))
order <- hclust(dist(as.matrix(temp)))$order

# Dont plot the T cell group from one patient
new_colnames <- c("Duodenum-B Cells", "Proximal-jejunum-B Cells", "Ileum-B Cells", "Colon-B Cells",
  "Duodenum-Plasma", "Proximal-jejunum-Plasma", "Ileum-Plasma", "Colon-Plasma",
  "Duodenum-CD4", "Proximal-jejunum-CD4", "Ileum-CD4", "Colon-CD4",
  "Duodenum-CD8", "Proximal-jejunum-CD8", "Ileum-CD8", "Colon-CD8",
  "Duodenum-NK", "Proximal-jejunum-NK", "Ileum-NK", "Colon-NK",
  "Duodenum-CyclingImmune", "Proximal-jejunum-CyclingImmune", "Ileum-CyclingImmune", "Colon-CyclingImmune",
  "Duodenum-Mast", "Proximal-jejunum-Mast", "Ileum-Mast", "Colon-Mast",
  "Duodenum-Mono_Macrophages", "Proximal-jejunum-Mono_Macrophages", "Ileum-Mono_Macrophages", "Colon-Mono_Macrophages",
  "DC",
  "ILC")#,
  #"T Cells")
   
test <- test[,new_colnames]
test <- .rowZscores(test)
ht1 <- .ArchRHeatmap(
    mat = test[order,],
    limits = c(-2,2),
    color = paletteContinuous(set = "blueYellow", n = 200), 
    clusterCols = FALSE, 
    clusterRows = FALSE,
    labelRows = TRUE,
    labelCols = TRUE,
    customColLabel = NULL,
    showRowDendrogram = FALSE,
    showColDendrogram = FALSE,
    draw = FALSE,
    #split = clustering$cluster, #scale = TRUE,
    name = paste0(nrow(test), " features\n")
)
plotPDF(ht1, name = paste0("Gene-Integration-TF-regulators-heatmap-", "all_sections", corr_savename, "corrcut", ".pdf"), width = 8, height = 15,  ArchRProj = proj_immune, addDOC = FALSE)

test_MM <- test_MM[,new_colnames]
test_MM <- .rowZscores(test_MM)
ht1 <- .ArchRHeatmap(
    mat = test_MM[order,],
    limits = c(-2,2),
    color = paletteContinuous(set = "solarExtra", n = 200), 
    clusterCols = FALSE, 
    clusterRows = FALSE,
    labelRows = TRUE,
    labelCols = TRUE,
    customColLabel = NULL,
    showRowDendrogram = FALSE,
    showColDendrogram = FALSE,
    draw = FALSE,
    #split = clustering$cluster, #scale = TRUE,
    name = paste0(nrow(test_MM), " features\n")
)
plotPDF(ht1, name = paste0("Motif-Matrix-TF-regulators-heatmap-", "all_sections", corr_savename, "corrcut", ".pdf"), width = 8, height = 15,  ArchRProj = proj_immune, addDOC = FALSE)


write.csv(test, paste0(parent_directory, "scATAC/projects/HuBMAP_immune_cells_multiome/Plots/Gene-Expression-Matrix-TF-regulators-heatmap-all_sections0p5corrcut.csv"))
write.csv(test_MM, paste0(parent_directory, "scATAC/projects/HuBMAP_immune_cells_multiome/Plots/Motif-Matrix-TF-regulators-heatmap-all_sections0p5corrcut.csv"))

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Stromal Cells
cellTypes <- getCellColData(proj_stromal)[,c("CellTypeRNA"), drop = FALSE]
cellLocations <- getCellColData(proj_stromal)[,c("Location"), drop = FALSE]
cellLocations[cellLocations$Location == "Transverse",] <- "Colon"
cellLocations[cellLocations$Location == "Sigmoid",] <- "Colon"
cellLocations[cellLocations$Location == "Descending",] <- "Colon"
cellLocations[cellLocations$Location == "Ascending",] <- "Colon"
cellLocations[cellLocations$Location == "Mid-jejunum",] <- "Proximal-jejunum"

cellTypes$CellTypeRNA <- paste0(cellLocations$Location, "-", cellTypes$CellTypeRNA)
cellTypes[cellTypes$CellTypeRNA == "Duodenum-Adipocytes",] <- "Adipocytes"
cellTypes[cellTypes$CellTypeRNA == "Proximal-jejunum-Adipocytes",] <- "Adipocytes"
cellTypes[cellTypes$CellTypeRNA == "Ileum-Adipocytes",] <- "Adipocytes"
cellTypes[cellTypes$CellTypeRNA == "Colon-Adipocytes",] <- "Adipocytes"

cellTypes[cellTypes$CellTypeRNA == "Duodenum-Neurons",] <- "Neurons"
cellTypes[cellTypes$CellTypeRNA == "Proximal-jejunum-Neurons",] <- "Neurons"
cellTypes[cellTypes$CellTypeRNA == "Ileum-Neurons",] <- "Neurons"
cellTypes[cellTypes$CellTypeRNA == "Colon-Neurons",] <- "Neurons"

cellTypes[cellTypes$CellTypeRNA == "Duodenum-ICC",] <- "ICC"
cellTypes[cellTypes$CellTypeRNA == "Proximal-jejunum-ICC",] <- "ICC"
cellTypes[cellTypes$CellTypeRNA == "Ileum-ICC",] <- "ICC"
cellTypes[cellTypes$CellTypeRNA == "Colon-ICC",] <- "ICC"

proj_stromal <- addCellColData(ArchRProj = proj_stromal, data = paste0(cellTypes$CellTypeRNA), cells = rownames(cellTypes), name = "CellTypeLocation", force = TRUE)

motifs_immune <- get_TF_regulators(proj_stromal = proj_stromal, saveName = "stromal", groupBy = "CellTypeLocation", corr = corr_cutoff)
motifs_immune <- motifs_immune[! duplicated(motifs_immune$GeneExpressionMatrix_name),]
motifs_immune <- motifs_immune[motifs_immune$GeneExpressionMatrix_name %ni% c("HNF1A-AS1", "RFX3-AS1", "MEIS1-AS2", "HNF4A-AS1", "MSC-AS1", "SOX2-OT", "TGIF2-RAB5IF", "CEBPB-AS1", "MEF2C-AS1", "NR2F1-AS1", "GATA2-AS1"),]

MM_epithelail <- getGroupSE(ArchRProj = proj_stromal, useMatrix = "VierstraMatrix", groupBy = "CellTypeLocation")
GIM_epithelail <- getGroupSE(ArchRProj = proj_stromal, useMatrix = "GeneExpressionMatrix", groupBy = "CellTypeLocation")

MM_epithelail_sub <- MM_epithelail[rowData(MM_epithelail)$name %in% motifs_immune$VierstraMatrix_name,]
GIM_epithelail_sub <- GIM_epithelail[rowData(GIM_epithelail)$name %in% motifs_immune$GeneExpressionMatrix_name,]

# Make a combined matrix with all cell types for plotting
test <- assays(GIM_epithelail_sub)$GeneExpressionMatrix
rownames(test) <- rowData(GIM_epithelail_sub)$name
colnames(test) <- colnames(GIM_epithelail_sub)
test <- test[motifs_immune$GeneExpressionMatrix_name,]

# get deviations, which is the second half of the matrix
test_MM <- assays(MM_epithelail_sub)$VierstraMatrix[(length(rowData(MM_epithelail_sub)$name)/2+1):length(rowData(MM_epithelail_sub)$name),]
rownames(test_MM) <- rowData(MM_epithelail_sub)$name[(length(rowData(MM_epithelail_sub)$name)/2+1):length(rowData(MM_epithelail_sub)$name)]
test_MM <- test_MM[motifs_immune$VierstraMatrix_name,]

# identify the order based on both the motif matrix and gene score matrix at the same time
temp <- cbind(.rowZscores(test), .rowZscores(test_MM))
order <- hclust(dist(as.matrix(temp)))$order

new_colnames <- c(
  "Duodenum-Myofibroblasts/SM 1", "Proximal-jejunum-Myofibroblasts/SM 1", "Ileum-Myofibroblasts/SM 1", "Colon-Myofibroblasts/SM 1",
  "Duodenum-Myofibroblasts/SM 2", "Proximal-jejunum-Myofibroblasts/SM 2", "Ileum-Myofibroblasts/SM 2", "Colon-Myofibroblasts/SM 2",
  "Duodenum-Myofibroblasts/SM 3", "Proximal-jejunum-Myofibroblasts/SM 3", "Ileum-Myofibroblasts/SM 3", "Colon-Myofibroblasts/SM 3",
  "Duodenum-Myofibroblasts/SM DES High", "Proximal-jejunum-Myofibroblasts/SM DES High", "Ileum-Myofibroblasts/SM DES High", "Colon-Myofibroblasts/SM DES High",
  "Duodenum-Pericytes", "Proximal-jejunum-Pericytes", "Ileum-Pericytes", "Colon-Pericytes",
  "Duodenum-Endothelial-Venules", "Proximal-jejunum-Endothelial-Venules", "Ileum-Endothelial-Venules", "Colon-Endothelial-Venules",
  "Duodenum-Endothelial-CD36+ Microvascular", "Proximal-jejunum-Endothelial-CD36+ Microvascular", "Ileum-Endothelial-CD36+ Microvascular", "Colon-Endothelial-CD36+ Microvascular",
  "Duodenum-Lymphatic Endothelial Cells", "Proximal-jejunum-Lymphatic Endothelial Cells", "Ileum-Lymphatic Endothelial Cells", "Colon-Lymphatic Endothelial Cells",
  "Duodenum-Crypt Fibroblasts 1 WNT2B+", "Proximal-jejunum-Crypt Fibroblasts 1 WNT2B+", "Ileum-Crypt Fibroblasts 1 WNT2B+", "Colon-Crypt Fibroblasts 1 WNT2B+",
  "Duodenum-Crypt Fibroblasts 2", "Proximal-jejunum-Crypt Fibroblasts 2", "Ileum-Crypt Fibroblasts 2", "Colon-Crypt Fibroblasts 2",
  "Duodenum-Crypt Fibroblasts 3 RSPO3+", "Proximal-jejunum-Crypt Fibroblasts 3 RSPO3+", "Ileum-Crypt Fibroblasts 3 RSPO3+", "Colon-Crypt Fibroblasts 3 RSPO3+",
  "Duodenum-Villus Fibroblasts WNT5B+", "Proximal-jejunum-Villus Fibroblasts WNT5B+", "Ileum-Villus Fibroblasts WNT5B+", "Colon-Villus Fibroblasts WNT5B+",
  "Duodenum-Glia", "Proximal-jejunum-Glia", "Ileum-Glia", "Colon-Glia",
  "Neurons",
  "Adipocytes",
  "ICC")

test <- test[,new_colnames]
test <- .rowZscores(test)
ht1 <- .ArchRHeatmap(
    mat = test[order,],
    limits = c(-2,2),
    color = paletteContinuous(set = "blueYellow", n = 200), 
    clusterCols = FALSE, 
    clusterRows = FALSE,
    labelRows = TRUE,
    labelCols = TRUE,
    customColLabel = NULL,
    showRowDendrogram = FALSE,
    showColDendrogram = FALSE,
    draw = FALSE,
    #split = clustering$cluster, #scale = TRUE,
    name = paste0(nrow(test), " features\n")
)
plotPDF(ht1, name = paste0("Gene-Integration-TF-regulators-heatmap-", "all_sections", corr_savename, "corrcut", ".pdf"), width = 10, height = 15,  ArchRProj = proj_stromal, addDOC = FALSE)

test_MM <- test_MM[,new_colnames]
test_MM <- .rowZscores(test_MM)
ht1 <- .ArchRHeatmap(
    mat = test_MM[order,],
    limits = c(-2,2),
    color = paletteContinuous(set = "solarExtra", n = 200), 
    clusterCols = FALSE, 
    clusterRows = FALSE,
    labelRows = TRUE,
    labelCols = TRUE,
    customColLabel = NULL,
    showRowDendrogram = FALSE,
    showColDendrogram = FALSE,
    draw = FALSE,
    #split = clustering$cluster, #scale = TRUE,
    name = paste0(nrow(test_MM), " features\n")
)
plotPDF(ht1, name = paste0("Motif-Matrix-TF-regulators-heatmap-", "all_sections", corr_savename, "corrcut", ".pdf"), width = 10, height = 15,  ArchRProj = proj_stromal, addDOC = FALSE)


write.csv(test, paste0(parent_directory, "scATAC/projects/HuBMAP_stromal_cells_multiome/Plots/Gene-Expression-Matrix-TF-regulators-heatmap-all_sections0p5corrcut.csv"))
write.csv(test_MM, paste0(parent_directory, "scATAC/projects/HuBMAP_stromal_cells_multiome/Plots/Motif-Matrix-TF-regulators-heatmap-all_sections0p5corrcut.csv"))

