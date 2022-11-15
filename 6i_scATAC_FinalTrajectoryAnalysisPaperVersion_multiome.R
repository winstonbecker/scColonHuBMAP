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

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Load previously defined archr projects
proj_ileum <- loadArchRProject(path = "HuBMAP_ileum_epithelial_cells_multiome")
proj_colon <- loadArchRProject(path = "HuBMAP_colon_epithelial_cells_multiome")
proj_duodenum <- loadArchRProject(path = "HuBMAP_duodenum_epithelial_cells_multiome")
proj_jejunum <- loadArchRProject(path = "HuBMAP_jejunum_epithelial_cells_multiome")
proj_epithelial <- loadArchRProject(path = "HuBMAP_epithelial_cells_multiome_final")

proj_duodenum <- addImputeWeights(proj_duodenum, reducedDims = paste("Harmony", "epithelial_duodenum", sep = ""))#, sampleCells = floor(nCells(proj_stromal)))
proj_jejunum <- addImputeWeights(proj_jejunum, reducedDims = paste("Harmony", "epithelial_jejunum", sep = ""))#, sampleCells = floor(nCells(proj_stromal)))
proj_ileum <- addImputeWeights(proj_ileum, reducedDims = paste("Harmony", "epithelial_ileum", sep = ""))#, sampleCells = floor(nCells(proj_stromal)))
proj_colon <- addImputeWeights(proj_colon, reducedDims = paste("Harmony", "epithelial_colon", sep = ""))#, sampleCells = floor(nCells(proj_stromal)))
proj_epithelial <- addImputeWeights(proj_epithelial, reducedDims = paste("Harmony", "hubmap_all", sep = ""))#, sampleCells = floor(nCells(proj_stromal)))

# Create peak trajectories
trajPMduodenum  <- getTrajectory(ArchRProj = proj_duodenum, name = "AbsorptiveU", useMatrix = "PeakMatrix", log2Norm = TRUE)
topPMduodenum <- plotTrajectoryHeatmap(trajPMduodenum,  varCutOff = 0.9,  returnMatrix = TRUE, scaleRows = FALSE, maxFeatures = 100000)
trajPMjejunum  <- getTrajectory(ArchRProj = proj_jejunum, name = "AbsorptiveU", useMatrix = "PeakMatrix", log2Norm = TRUE)
topPMjejunum <- plotTrajectoryHeatmap(trajPMjejunum,  varCutOff = 0.9,  returnMatrix = TRUE, scaleRows = FALSE, maxFeatures = 100000)
trajPMileum  <- getTrajectory(ArchRProj = proj_ileum, name = "AbsorptiveU", useMatrix = "PeakMatrix", log2Norm = TRUE)
topPMileum <- plotTrajectoryHeatmap(trajPMileum,  varCutOff = 0.9,  returnMatrix = TRUE, scaleRows = FALSE, maxFeatures = 100000)
trajPMcolon  <- getTrajectory(ArchRProj = proj_colon, name = "AbsorptiveU", useMatrix = "PeakMatrix", log2Norm = TRUE)
topPMcolon <- plotTrajectoryHeatmap(trajPMcolon,  varCutOff = 0.9,  returnMatrix = TRUE, scaleRows = FALSE, maxFeatures = 100000)

# Create gene integration matrix trajectories
trajGIMduodenum  <- getTrajectory(ArchRProj = proj_duodenum, name = "AbsorptiveU", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
topGIMduodenum <- plotTrajectoryHeatmap(trajGIMduodenum,  varCutOff = 0.9,  returnMatrix = TRUE, scaleRows = FALSE, maxFeatures = 100000)
trajGIMjejunum  <- getTrajectory(ArchRProj = proj_jejunum, name = "AbsorptiveU", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
topGIMjejunum <- plotTrajectoryHeatmap(trajGIMjejunum,  varCutOff = 0.9,  returnMatrix = TRUE, scaleRows = FALSE, maxFeatures = 100000)
trajGIMileum  <- getTrajectory(ArchRProj = proj_ileum, name = "AbsorptiveU", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
topGIMileum <- plotTrajectoryHeatmap(trajGIMileum,  varCutOff = 0.9,  returnMatrix = TRUE, scaleRows = FALSE, maxFeatures = 100000)
trajGIMcolon  <- getTrajectory(ArchRProj = proj_colon, name = "AbsorptiveU", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
topGIMcolon <- plotTrajectoryHeatmap(trajGIMcolon,  varCutOff = 0.9,  returnMatrix = TRUE, scaleRows = FALSE, maxFeatures = 100000)

# Create motif matrix trajectories
trajMMduodenum  <- getTrajectory(ArchRProj = proj_duodenum, name = "AbsorptiveU", useMatrix = "VierstraMatrix", log2Norm = TRUE)
trajMMjejunum  <- getTrajectory(ArchRProj = proj_jejunum, name = "AbsorptiveU", useMatrix = "VierstraMatrix", log2Norm = TRUE)
trajMMileum  <- getTrajectory(ArchRProj = proj_ileum, name = "AbsorptiveU", useMatrix = "VierstraMatrix", log2Norm = TRUE)
trajMMcolon  <- getTrajectory(ArchRProj = proj_colon, name = "AbsorptiveU", useMatrix = "VierstraMatrix", log2Norm = TRUE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Peak matrix analysis
# Add numbers to trajectory columns to maintain order
colnames(trajPMduodenum) <- paste0("1.", colnames(trajPMduodenum))
colnames(trajPMjejunum) <- paste0("2.", colnames(trajPMjejunum))
colnames(trajPMileum) <- paste0("3.", colnames(trajPMileum))
colnames(trajPMcolon) <- paste0("4.", colnames(trajPMcolon))

fullTrajPM <- cbind((trajPMduodenum), (trajPMjejunum), (trajPMileum), (trajPMcolon))
fullTrajPM <- fullTrajPM[unique(c(rownames(topPMduodenum), rownames(topPMcolon), rownames(topPMjejunum), rownames(topPMileum))),]
fullTrajPMtop <- fullTrajPM[(rowMax(assays(fullTrajPM)$mat)-rowMin(assays(fullTrajPM)$mat))>0.2,] # only keep things with an absolute difference of 0.2
test <- .rowZscores(assays(fullTrajPMtop)$mat)
# Cluster and plot
set.seed(1) # set the seed
clustering <- kmeans(test, 7, iter.max = 500)
ht1 <- .ArchRHeatmap(
    mat = test,
    limits = c(-2,2),
    color = paletteContinuous(set = "solarExtra", n = 200), 
    clusterCols = FALSE, 
    clusterRows = TRUE,
    labelRows = FALSE,
    labelCols = FALSE,
    customColLabel = NULL,
    showRowDendrogram = FALSE,
    showColDendrogram = FALSE,
    draw = FALSE,
    split = clustering$cluster, #scale = TRUE,
    name = paste0(nrow(test), " features\n")
)
plotPDF(ht1, name = paste0("full_clustered_peak_trajectory-r.pdf"), width = 8, height = 12, ArchRProj = proj_epithelial, addDOC = FALSE)

temp <-  clustering$cluster
temp[temp == 5] <- "bad"
temp[temp == 1] <- "bad"
temp[temp == 4] <- "LateColon"
temp[temp == 7] <- "Late"
temp[temp == 3] <- "LateSmallBowel"
temp[temp == 2] <- "Early"
temp[temp == 6] <- "EarlyColon"
test2<-as.data.frame(test)
test2$cluster <- temp
test2 <- test2[test2$cluster != "bad",]
write.table(format(test2, digits = 2,nsmall = 3,scientific=F), paste0(parent_directory, "scATAC/projects/HuBMAP_epithelial_cells_multiome_final/Plots/Figure6B_peak_trajectories.tsv"),sep = "\t", quote = FALSE)
write.table(test2[,"cluster",drop = FALSE], paste0(parent_directory, "scATAC/projects/HuBMAP_epithelial_cells_multiome_final/Plots/Figure6B_peak_trajectories_cluster_ONLY.tsv"),sep = "\t", quote = FALSE)


# Identify motifs enriched in the clusters of peaks
# create dummy marker test 
markerTest <- getMarkerFeatures(
  ArchRProj = proj_colon, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeRNA",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Stem",
  bgdGroups = "Enterocytes"
)
# iterate through clusters to find enriched motifs
nclust <- 7
for (i in 1:nclust){
    clusters <- clustering$cluster
    clusters <- clusters[clusters == i]
    clusterMarkerTest <- markerTest
    # now just set the peaks in the cluster as the only ones that will meet the significance cutoff
    assays(clusterMarkerTest)$Log2FC[,] <- 0
    assays(clusterMarkerTest)$Log2FC[paste0(rowData(clusterMarkerTest)$seqnames, ":", rowData(clusterMarkerTest)$start, "_", rowData(clusterMarkerTest)$end) %in% names(clusters),] <- 10000
    motifsUp <- peakAnnoEnrichment(
      seMarker = clusterMarkerTest,
      ArchRProj = proj_colon,
      peakAnnotation = "Vierstra",
      cutOff = "Log2FC >9999"
    )

    df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
    df <- df[order(df$mlog10Padj, decreasing = TRUE),]
    df$rank <- seq_len(nrow(df))

    new_additionUp <- assays(motifsUp)$mlog10Padj
    colnames(new_additionUp) = paste0("C", i)
    if (i==1){
      save_motifs_up <- new_additionUp
    }
    if (i>1){
      save_motifs_up <- cbind(new_additionUp, save_motifs_up)
    }
}
df <- DataFrame(save_motifs_up)
motif_names <- names(apply(save_motifs_up, 1, max)[order(apply(save_motifs_up, 1, max), decreasing = TRUE)])
motif_names_df <- t(data.frame(strsplit(motif_names, ":")))
motif_names <- motif_names[!duplicated(motif_names_df[,3])]
#non-redundant df
motif_names_df <- motif_names_df[!duplicated(motif_names_df[,3]),]
df <- DataFrame(save_motifs_up[motif_names[order(as.numeric(motif_names_df[,3]))],])
df <- df[apply(df, 1, function(x) !all(x==0)),]
df <- df[,paste0("C", c(2:4,6,7))]
df <- df[apply(df, 1, function(x) sum(x>75)>0),]

paletteLength <- 256
myBreaks <- c(seq(0, 50, length.out=ceiling(paletteLength/2) + 1), 
            seq(50.1, 200, length.out=floor(paletteLength/2)))
set.seed(100)
p <- pheatmap::pheatmap(
  mat = as.matrix(df), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", breaks = myBreaks, cluster_cols = TRUE, cluster_rows = FALSE
)
plotPDF(p, name = paste0("test_Motif-Enrichment_Cluster-Trajectory"), width = 20, height = 30,  ArchRProj = proj_epithelial, addDOC = FALSE)



############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Gene Expression Matrix Analysis
markerGenes <- c("CA2", "SI", # absorptive
    "SOX9", #progenitor
    "SMOC2", "RGMB", "LGR5", "ASCL2", #stem
    "HNF4A", "HNF4G", "GPX2",
    "TMPRSS15", "ETV6", 
    "CLCA4",
    "PLB1", "CUBN", # Ileum
    "CEACAM5", "SLC26A3",
    "SLC4A7", "APOB",
    "REG3A", "SCTR", "MTTP", "APOA4", "APOC3", "MME", "ROR2", "PTGDR", "SCNN1B", "P3H2", "MS4A12") #LIPID

colnames(trajGIMduodenum) <- paste0("1.", colnames(trajGIMduodenum))
colnames(trajGIMjejunum) <- paste0("2.", colnames(trajGIMjejunum))
colnames(trajGIMileum) <- paste0("3.", colnames(trajGIMileum))
colnames(trajGIMcolon) <- paste0("4.", colnames(trajGIMcolon))

fullTrajGIM <- cbind((trajGIMduodenum), (trajGIMjejunum), (trajGIMileum), (trajGIMcolon))
fullTrajGIM <- fullTrajGIM[unique(c(rownames(topGIMduodenum), rownames(topGIMcolon), rownames(topGIMjejunum), rownames(topGIMileum))),]
fullTrajGIMtop <- fullTrajGIM[(rowMax(assays(fullTrajGIM)$mat)-rowMin(assays(fullTrajGIM)$mat))>0.5,]
test <- .rowZscores(assays(fullTrajGIMtop)$mat)
set.seed(1)
clustering <- kmeans(test, 7, iter.max = 500)
idxLabel <- unique(unlist(lapply(paste0(":",markerGenes), function(x) rownames(test)[grepl(x, rownames(test))])))

ht1 <- .ArchRHeatmap(
    mat = test,
    limits = c(-1.75,1.75),
    color = paletteContinuous(set = "blueYellow", n = 200), 
    clusterCols = FALSE, 
    clusterRows = TRUE,
    labelRows = FALSE,
    labelCols = FALSE,
    customColLabel = NULL,
    showRowDendrogram = FALSE,
    showColDendrogram = FALSE,
    draw = FALSE,
    split = clustering$cluster, #scale = TRUE,
    customRowLabel = match(idxLabel, 
                rownames(test)),
    name = paste0(nrow(test), " features\n")
)
plotPDF(ht1, name = paste0("full_kmeans_clustered_rowmax_ordered_GIM_trajectory-r.pdf"), width = 8, height = 12, ArchRProj = proj_epithelial, addDOC = FALSE)

saveRDS(clustering, "clustering_genes_traj.rds")

temp <-  clustering$cluster
temp[temp == 1] <- "EarlyDuodenumJejunum"
temp[temp == 2] <- "Ileum"
temp[temp == 3] <- "LateDuodenumJejunum"
temp[temp == 4] <- "EarlyColon"
temp[temp == 5] <- "Early"
temp[temp == 6] <- "LateColon"
temp[temp == 7] <- "LateSmallIntestine"
test2<-as.data.frame(test)
test2$cluster <- temp
write.table(format(test2, digits = 2,nsmall = 3,scientific=F), paste0(parent_directory, "scATAC/projects/HuBMAP_epithelial_cells_multiome_final/Plots/Figure6C_gene_trajectories.tsv"),sep = "\t", quote = FALSE)
write.table(test2[,"cluster",drop = FALSE], paste0(parent_directory, "scATAC/projects/HuBMAP_epithelial_cells_multiome_final/Plots/Figure6C_gene_trajectories_cluster_ONLY.tsv"),sep = "\t", quote = FALSE)

fullTrajGIMtop <- fullTrajGIM[(rowMax(assays(fullTrajGIM)$mat)-rowMin(assays(fullTrajGIM)$mat))>4,]
test <- .rowZscores(assays(fullTrajGIMtop)$mat)
set.seed(1)
clustering <- kmeans(test, 7, iter.max = 500)
idxLabel <- unique(unlist(lapply(rownames(fullTrajGIMtop), function(x) rownames(test)[grepl(x, rownames(test))])))

ht1 <- .ArchRHeatmap(
    mat = test,
    limits = c(-2,2),
    color = paletteContinuous(set = "blueYellow", n = 200), 
    clusterCols = FALSE, 
    clusterRows = TRUE,
    labelRows = FALSE,
    labelCols = FALSE,
    customColLabel = NULL,
    showRowDendrogram = FALSE,
    showColDendrogram = FALSE,
    draw = FALSE,
    split = clustering$cluster, #scale = TRUE,
    customRowLabel = match(idxLabel, 
                rownames(test)),
    name = paste0(nrow(test), " features\n")
)
plotPDF(ht1, name = paste0("full_kmeans_clustered_rowmax_ordered_GIM_trajectory_l2fc_gt_4.pdf"), width = 8, height = 12, ArchRProj = proj_epithelial, addDOC = FALSE)



##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# correlate GIM and MM along trajectories
find_correlated_trajectories <- function(proj, trajGIMduodenum, trajMMduodenum){
    # Look for ones that are correlated along the trajectory
    corGIM_MM_d <- correlateTrajectories(trajGIMduodenum, trajMMduodenum)
    corGIM_MM_d[[1]]
    trajGIMduodenum2 <- trajGIMduodenum[corGIM_MM_d[[1]]$name1, ]
    trajMMduodenum2 <- trajMMduodenum[corGIM_MM_d[[1]]$name2, ]

    # same gene and same group
    keep <- !(duplicated(rownames(trajGIMduodenum2)) & duplicated(t(data.frame(strsplit(rownames(trajMMduodenum2), ":"))[4,])))
    trajGIMduodenum2 <- trajGIMduodenum2[keep,]
    trajMMduodenum2 <- trajMMduodenum2[keep,]

    # define ordering
    trajCombined <- trajGIMduodenum2
    p0 <- t(apply(assay(trajGIMduodenum2), 1, scale))
    p1 <- t(apply(assay(trajMMduodenum2), 1, scale))
    rownames(p1) <- rownames(p0)
    p0 <- p0+p1
    dimnames(p0) <- dimnames(assay(trajCombined))
    assay(trajCombined) <- p0

    combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
    rowOrder <- match(rownames(combinedMat), rownames(trajGIMduodenum2))

    ht1 <- plotTrajectoryHeatmap(trajGIMduodenum2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelRows = TRUE)#labelTop = 5000)
    ht2 <- plotTrajectoryHeatmap(trajMMduodenum2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelRows = TRUE, scaleRows = FALSE)#, labelTop = 5000)
    plotPDF(ht1, ht2, name = "Correlated_trajectory.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 12)

    return(list(trajGIMduodenum2[rowOrder,], trajMMduodenum2[rowOrder,]))
}

test <- find_correlated_trajectories(proj_duodenum, trajGIMduodenum, trajMMduodenum)
trajGIMduodenum2 <- test[[1]]
trajMMduodenum2 <- test[[2]]

test <- find_correlated_trajectories(proj_jejunum, trajGIMjejunum, trajMMjejunum)
trajGIMjejunum2 <- test[[1]]
trajMMjejunum2 <- test[[2]]

test <- find_correlated_trajectories(proj_ileum, trajGIMileum, trajMMileum)
trajGIMileum2 <- test[[1]]
trajMMileum2 <- test[[2]]

test <- find_correlated_trajectories(proj_colon, trajGIMcolon, trajMMcolon)
trajGIMcolon2 <- test[[1]]
trajMMcolon2 <- test[[2]]


motifs <- unique(c(rownames(trajGIMcolon2), rownames(trajGIMileum2), rownames(trajGIMjejunum2), rownames(trajGIMduodenum2)))
subscript = "hubmap_all"
proj_stromal <- proj_epithelial
proj_stromal <- addImputeWeights(proj_stromal, reducedDims = paste("Harmony", subscript, sep = ""))#, sampleCells = floor(nCells(proj_stromal)))

p <- plotEmbedding(
    ArchRProj = proj_stromal, 
    colorBy = "GeneExpressionMatrix", 
    name = motifs,
    pal = paletteContinuous("blueYellow"),
    embedding = paste("UMAP", subscript, sep = ""),
    imputeWeights = getImputeWeights(proj_stromal)
)
plotPDF(p, name = paste0("test-Gene-Integration-Correlated_Along_Trajectory-", subscript, ".pdf"), width = 5, height = 5,  ArchRProj = proj_stromal, addDOC = FALSE)


AbsorptiveU_duodenum <- getCellColData(proj_duodenum)[,"AbsorptiveU", drop = FALSE]
AbsorptiveU_duodenum <- AbsorptiveU_duodenum[!is.na(AbsorptiveU_duodenum),,drop = FALSE]
AbsorptiveU_jejunum <- getCellColData(proj_jejunum)[,"AbsorptiveU", drop = FALSE]
AbsorptiveU_jejunum <- AbsorptiveU_jejunum[!is.na(AbsorptiveU_jejunum),,drop = FALSE]
AbsorptiveU_ileum <- getCellColData(proj_ileum)[,"AbsorptiveU", drop = FALSE]
AbsorptiveU_ileum <- AbsorptiveU_ileum[!is.na(AbsorptiveU_ileum),,drop = FALSE]
AbsorptiveU_colon <- getCellColData(proj_colon)[,"AbsorptiveU", drop = FALSE]
AbsorptiveU_colon <- AbsorptiveU_colon[!is.na(AbsorptiveU_colon),,drop = FALSE]

proj_stromal <- addCellColData(ArchRProj = proj_stromal,data = AbsorptiveU_duodenum$AbsorptiveU,name = "AbsorptiveU_duodenum",cells = rownames(AbsorptiveU_duodenum),force = FALSE)
proj_stromal <- addCellColData(ArchRProj = proj_stromal,data = AbsorptiveU_jejunum$AbsorptiveU,name = "AbsorptiveU_jejunum",cells = rownames(AbsorptiveU_jejunum),force = FALSE)
proj_stromal <- addCellColData(ArchRProj = proj_stromal,data = AbsorptiveU_ileum$AbsorptiveU,name = "AbsorptiveU_ileum",cells = rownames(AbsorptiveU_ileum),force = FALSE)
proj_stromal <- addCellColData(ArchRProj = proj_stromal,data = AbsorptiveU_colon$AbsorptiveU,name = "AbsorptiveU_colon",cells = rownames(AbsorptiveU_colon),force = FALSE)


trajGIMduodenum_from_full_proj  <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU_duodenum", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
trajMMduodenum_from_full_proj  <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU_duodenum", useMatrix = "VierstraMatrix", log2Norm = FALSE)
trajGIMjejunum_from_full_proj  <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU_jejunum", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
trajMMjejunum_from_full_proj  <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU_jejunum", useMatrix = "VierstraMatrix", log2Norm = FALSE)
trajGIMileum_from_full_proj  <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU_ileum", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
trajMMileum_from_full_proj  <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU_ileum", useMatrix = "VierstraMatrix", log2Norm = FALSE)
trajGIMcolon_from_full_proj  <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU_colon", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
trajMMcolon_from_full_proj  <- getTrajectory(ArchRProj = proj_stromal, name = "AbsorptiveU_colon", useMatrix = "VierstraMatrix", log2Norm = FALSE)

colnames(trajGIMduodenum_from_full_proj) <- paste0("1.", colnames(trajGIMduodenum_from_full_proj))
colnames(trajGIMjejunum_from_full_proj) <- paste0("2.", colnames(trajGIMjejunum_from_full_proj))
colnames(trajGIMileum_from_full_proj) <- paste0("3.", colnames(trajGIMileum_from_full_proj))
colnames(trajGIMcolon_from_full_proj) <- paste0("4.", colnames(trajGIMcolon_from_full_proj))


motifs <- unique(c(rownames(trajMMcolon2), rownames(trajMMileum2), rownames(trajMMjejunum2), rownames(trajMMduodenum2)))
genes <- unique(c(rownames(trajGIMcolon2), rownames(trajGIMileum2), rownames(trajGIMjejunum2), rownames(trajGIMduodenum2)))


full_mat <- cbind(assays(trajGIMduodenum_from_full_proj[genes,6:95])$smoothMat,
      assays(trajGIMjejunum_from_full_proj[genes,6:95])$smoothMat,
      assays(trajGIMileum_from_full_proj[genes,6:95])$smoothMat,
      assays(trajGIMcolon_from_full_proj[genes,6:95])$smoothMat)


#order <- names(rowSums(full_mat*c(1:100,1:100,1:100,1:100))[order(rowSums(.rowZscores(full_mat)%*%c(1:100,1:100,1:100,1:100)))])
order <- rownames(full_mat)[order(.rowZscores(full_mat)%*%(c((6:95)**1.5,(6:95)**1.5,(6:95)**1.5,(6:95)**1.5)))]
library(stringr)
motif_order <- c()
for (gene in order){
  motif_order <- c(motif_order, motifs[grepl(paste0(str_split(gene, ":")[[1]][2],"_"), motifs)])
}
# normalize for number of cells (otherwise seems to just be a sum)
full_mat_motifs <- cbind(assays(trajMMduodenum_from_full_proj[motif_order,6:95])$smoothMat/length(AbsorptiveU_duodenum$AbsorptiveU),
      assays(trajMMjejunum_from_full_proj[motif_order,6:95])$smoothMat/length(AbsorptiveU_jejunum$AbsorptiveU),
      assays(trajMMileum_from_full_proj[motif_order,6:95])$smoothMat/length(AbsorptiveU_ileum$AbsorptiveU),
      assays(trajMMcolon_from_full_proj[motif_order,6:95])$smoothMat/length(AbsorptiveU_colon$AbsorptiveU))

order <- order[order %ni% c("chr12:HNF1A-AS1", "chr9:RFX3-AS1", "chr20:HNF4A-AS1", "chr17:SOX9-AS1", "chr3:SOX2-OT", "chr1:ELF3-AS1")]

test <- .rowZscores(full_mat[order,])
ht1 <- .ArchRHeatmap(
    mat = test,
    limits = c(-2,2),
    color = paletteContinuous(set = "blueYellow", n = 200), 
    clusterCols = FALSE, 
    clusterRows = FALSE,
    labelRows = TRUE,
    labelCols = FALSE,
    customColLabel = NULL,
    showRowDendrogram = FALSE,
    showColDendrogram = FALSE,
    draw = FALSE,
    name = paste0(nrow(test), " features\n")
)
plotPDF(ht1, name = paste0("integrated-GIM-matrix-AS1-removed-reproduce.pdf"), width = 8, height = 12, ArchRProj = proj_epithelial, addDOC = FALSE)

write.table(test, paste0(parent_directory, "scATAC/projects/HuBMAP_epithelial_cells_multiome_final/Plots/Figure6D_correlated_trajectories.tsv"),sep = "\t", quote = FALSE)

##############################################################################################################
#............................................................................................................................#
##############################################################################################################################
# plot some examples
proj_duodenum <- addImputeWeights(proj_duodenum, reducedDims = paste("Harmony", "epithelial_duodenum", sep = ""))#, sampleCells = floor(nCells(proj_stromal)))
proj_jejunum <- addImputeWeights(proj_jejunum, reducedDims = paste("Harmony", "epithelial_jejunum", sep = ""))#, sampleCells = floor(nCells(proj_stromal)))
proj_ileum <- addImputeWeights(proj_ileum, reducedDims = paste("Harmony", "epithelial_ileum", sep = ""))#, sampleCells = floor(nCells(proj_stromal)))
proj_colon <- addImputeWeights(proj_colon, reducedDims = paste("Harmony", "epithelial_colon", sep = ""))#, sampleCells = floor(nCells(proj_stromal)))

p1 <- plotTrajectory(proj_duodenum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_duodenum"), colorBy = "GeneExpressionMatrix", name = "ETV6", continuousSet = "horizonExtra")
plotPDF(p1, name = "ETV6-Plot-duodenum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_jejunum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_jejunum"), colorBy = "GeneExpressionMatrix", name = "ETV6", continuousSet = "horizonExtra")
plotPDF(p1, name = "ETV6-Plot-jejunum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_ileum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_ileum"), colorBy = "GeneExpressionMatrix", name = "ETV6", continuousSet = "horizonExtra")
plotPDF(p1, name = "ETV6-Plot-ileum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_colon, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_colon"), colorBy = "GeneExpressionMatrix", name = "ETV6", continuousSet = "horizonExtra")
plotPDF(p1, name = "ETV6-Plot-colon_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)

p1 <- plotTrajectory(proj_duodenum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_duodenum"), colorBy = "GeneExpressionMatrix", name = "TMPRSS15", continuousSet = "horizonExtra")
plotPDF(p1, name = "TMPRSS15-Plot-duodenum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_jejunum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_jejunum"), colorBy = "GeneExpressionMatrix", name = "TMPRSS15", continuousSet = "horizonExtra")
plotPDF(p1, name = "TMPRSS15-Plot-jejunum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_ileum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_ileum"), colorBy = "GeneExpressionMatrix", name = "TMPRSS15", continuousSet = "horizonExtra")
plotPDF(p1, name = "TMPRSS15-Plot-ileum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_colon, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_colon"), colorBy = "GeneExpressionMatrix", name = "TMPRSS15", continuousSet = "horizonExtra")
plotPDF(p1, name = "TMPRSS15-Plot-colon_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)

p1 <- plotTrajectory(proj_duodenum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_duodenum"), colorBy = "GeneExpressionMatrix", name = "LGR5", continuousSet = "horizonExtra")
plotPDF(p1, name = "LGR5-Plot-duodenum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_jejunum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_jejunum"), colorBy = "GeneExpressionMatrix", name = "LGR5", continuousSet = "horizonExtra")
plotPDF(p1, name = "LGR5-Plot-jejunum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_ileum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_ileum"), colorBy = "GeneExpressionMatrix", name = "LGR5", continuousSet = "horizonExtra")
plotPDF(p1, name = "LGR5-Plot-ileum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_colon, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_colon"), colorBy = "GeneExpressionMatrix", name = "LGR5", continuousSet = "horizonExtra")
plotPDF(p1, name = "LGR5-Plot-colon_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)

p1 <- plotTrajectory(proj_duodenum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_duodenum"), colorBy = "GeneExpressionMatrix", name = "SCNN1B", continuousSet = "horizonExtra")
plotPDF(p1, name = "SCNN1B-Plot-duodenum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_jejunum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_jejunum"), colorBy = "GeneExpressionMatrix", name = "SCNN1B", continuousSet = "horizonExtra")
plotPDF(p1, name = "SCNN1B-Plot-jejunum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_ileum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_ileum"), colorBy = "GeneExpressionMatrix", name = "SCNN1B", continuousSet = "horizonExtra")
plotPDF(p1, name = "SCNN1B-Plot-ileum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_colon, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_colon"), colorBy = "GeneExpressionMatrix", name = "SCNN1B", continuousSet = "horizonExtra")
plotPDF(p1, name = "SCNN1B-Plot-colon_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)

p1 <- plotTrajectory(proj_duodenum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_duodenum"), colorBy = "GeneExpressionMatrix", name = "MTTP", continuousSet = "horizonExtra")
plotPDF(p1, name = "MTTP-Plot-duodenum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_jejunum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_jejunum"), colorBy = "GeneExpressionMatrix", name = "MTTP", continuousSet = "horizonExtra")
plotPDF(p1, name = "MTTP-Plot-jejunum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_ileum, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_ileum"), colorBy = "GeneExpressionMatrix", name = "MTTP", continuousSet = "horizonExtra")
plotPDF(p1, name = "MTTP-Plot-ileum_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
p1 <- plotTrajectory(proj_colon, trajectory = "AbsorptiveU", embedding = paste0("UMAPHarmonyepithelial_colon"), colorBy = "GeneExpressionMatrix", name = "MTTP", continuousSet = "horizonExtra")
plotPDF(p1, name = "MTTP-Plot-colon_abs_diff-Genes-Traj-UMAP.pdf", ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)

##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# look for correlated peak changes with enhancer changes

get_smoothed_pseudotime <- function(proj_duodenum, proj_jejunum, proj_ileum, proj_colon, traj, gene, correlation = 0.6){
    p2g <- getPeak2GeneLinks(ArchRProj = proj_duodenum, corCutOff = correlation, resolution = 1,returnLoops = FALSE)
    correlated_peaks <- metadata(p2g)[[1]][p2g[p2g$idxRNA == which(metadata(p2g)[[2]]$name == gene),"idxATAC"]]
    correlated_peaks <- paste0(seqnames(correlated_peaks), ":", start(ranges(correlated_peaks)), "_", end(ranges(correlated_peaks)))

    p2g <- getPeak2GeneLinks(ArchRProj = proj_jejunum, corCutOff = correlation, resolution = 1,returnLoops = FALSE)
    correlated_peaks_new <- metadata(p2g)[[1]][p2g[p2g$idxRNA == which(metadata(p2g)[[2]]$name == gene),"idxATAC"]]
    correlated_peaks <- c(correlated_peaks, paste0(seqnames(correlated_peaks_new), ":", start(ranges(correlated_peaks_new)), "_", end(ranges(correlated_peaks_new))))

    p2g <- getPeak2GeneLinks(ArchRProj = proj_ileum, corCutOff = correlation, resolution = 1,returnLoops = FALSE)
    correlated_peaks_new <- metadata(p2g)[[1]][p2g[p2g$idxRNA == which(metadata(p2g)[[2]]$name == gene),"idxATAC"]]
    correlated_peaks <- c(correlated_peaks, paste0(seqnames(correlated_peaks_new), ":", start(ranges(correlated_peaks_new)), "_", end(ranges(correlated_peaks_new))))

    #correlated_peaks <- c()
    p2g <- getPeak2GeneLinks(ArchRProj = proj_colon, corCutOff = correlation, resolution = 1,returnLoops = FALSE)
    correlated_peaks_new <- metadata(p2g)[[1]][p2g[p2g$idxRNA == which(metadata(p2g)[[2]]$name == gene),"idxATAC"]]
    correlated_peaks <- unique(c(correlated_peaks, paste0(seqnames(correlated_peaks_new), ":", start(ranges(correlated_peaks_new)), "_", end(ranges(correlated_peaks_new)))))

    correlated_peaks <- correlated_peaks[correlated_peaks != ":_"]

    new <- assays(traj[correlated_peaks,])$smoothMat
    return(new)
}

plot_peak_pseudotime <- function(new, linked_row_maxes, saveName, gene, correlation = 0.6, ymax = 0.25){
    rownames_new <- rownames(new)
    col1 <- c()
    col2 <- c()
    col3 <- c()
    for (i in rownames(new)){
        # col1 <- c(col1, rep(i,100))
        # col2 <- c(col2, 0:99+0.5)
        # col3 <- c(col3, new[i,])
        col1 <- c(col1, rep(i,90))
        col2 <- c(col2, 5:94+0.5)
        col3 <- c(col3, new[i,6:95])
    }
    df <- data.frame(peak = col1, pseudotime = col2, access = col3)
    p1 <- ggplot(data=df, aes(x=pseudotime, y=access, group=peak)) + geom_line(aes(color=peak)) + geom_point(aes(color=peak))+scale_color_manual(values =paste0(ArchRPalettes$stallion)) + theme_ArchR() + theme(aspect.ratio=5/2) + scale_y_continuous(limits = c(0, ymax))
    plotPDF(p1, name = paste0(saveName, "-linked-peaks-", gene, "-vs-pseudotime-set-ymax.pdf"), ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)

    new <- diag(1/linked_row_maxes) %*% new
    rownames(new) <- rownames_new
    col1 <- c()
    col2 <- c()
    col3 <- c()
    for (i in rownames(new)){
        # col1 <- c(col1, rep(i,100))
        # col2 <- c(col2, 0:99+0.5)
        # col3 <- c(col3, new[i,])
        col1 <- c(col1, rep(i,90))
        col2 <- c(col2, 5:94+0.5)
        col3 <- c(col3, new[i,6:95])
    }
    df <- data.frame(peak = col1, pseudotime = col2, access = col3)
    p1 <- ggplot(data=df, aes(x=pseudotime, y=access, group=peak)) + geom_line(aes(color=peak)) + geom_point(aes(color=peak))+scale_color_manual(values =paste0(ArchRPalettes$stallion)) + theme_ArchR() + theme(aspect.ratio=5/2) + scale_y_continuous(limits = c(0, 1))
    plotPDF(p1, name = paste0(saveName, "-linked-peaks-", gene, "-vs-pseudotime-normalized-to-1.pdf"), ArchRProj = proj_epithelial, addDOC = FALSE, width = 5, height = 5)
    write.table(df, paste0(parent_directory, "scATAC/projects/HuBMAP_epithelial_cells_multiome_final/Plots/", saveName, "-linked-peaks-", gene, "-vs-pseudotime-normalized-to-1.tsv"),sep = "\t", quote = FALSE)
}

peaks <- getPeakSet(ArchRProj = proj_duodenum)

gene_list <- c("ETV6") 
correlation_list <- c(0.4, 0.4, 0.4, 0.4, 0.4)

gene_list <- c("TMPRSS15") 
correlation_list <- c(0.55, 0.55, 0.55, 0.55, 0.55)

for (j in 1:length(gene_list)){
  correlation <- correlation_list[j]
  gene <- gene_list[j]

  duod_linked <- get_smoothed_pseudotime(proj_duodenum, proj_jejunum, proj_ileum, proj_colon, traj = trajPMduodenum, gene = gene, correlation = correlation)
  jejunum_linked <- get_smoothed_pseudotime(proj_duodenum, proj_jejunum, proj_ileum, proj_colon, traj = trajPMjejunum, gene = gene, correlation = correlation)
  ileum_linked <- get_smoothed_pseudotime(proj_duodenum, proj_jejunum, proj_ileum, proj_colon, traj = trajPMileum, gene = gene, correlation = correlation)
  colon_linked <- get_smoothed_pseudotime(proj_duodenum, proj_jejunum, proj_ileum, proj_colon, traj = trajPMcolon, gene = gene, correlation = correlation)
  linked_row_maxes <- rowMax(cbind(duod_linked, jejunum_linked, ileum_linked, colon_linked))

  print(peaks[paste0(seqnames(peaks), ":", start(ranges(peaks)), "_", end(ranges(peaks))) %in% rownames(duod_linked),])

  # this section is just to add it it is intronic, promoter, etc.
  peak_sub <- data.frame(peaks[paste0(seqnames(peaks), ":", start(ranges(peaks)), "_", end(ranges(peaks))) %in% rownames(duod_linked),])
  rownames(peak_sub) <- paste0(peak_sub$seqnames, ":", peak_sub$start, "_", peak_sub$end)
  peak_sub$seq_peak_type <- paste0(peak_sub$seqnames, ":", peak_sub$start, "_", peak_sub$end, "_", peak_sub$peakType)
  rownames(duod_linked) <- peak_sub[rownames(duod_linked), "seq_peak_type"]
  rownames(jejunum_linked) <- peak_sub[rownames(jejunum_linked), "seq_peak_type"]
  rownames(ileum_linked) <- peak_sub[rownames(ileum_linked), "seq_peak_type"]
  rownames(colon_linked) <- peak_sub[rownames(colon_linked), "seq_peak_type"]

  plot_peak_pseudotime(duod_linked, linked_row_maxes, saveName = "duodenum", gene = gene, correlation = correlation)
  plot_peak_pseudotime(jejunum_linked, linked_row_maxes, saveName = "jejunum", gene = gene, correlation = correlation)
  plot_peak_pseudotime(ileum_linked, linked_row_maxes, saveName = "ileum", gene = gene, correlation = correlation)
  plot_peak_pseudotime(colon_linked, linked_row_maxes, saveName = "colon", gene = gene, correlation = correlation)
}





##############################################################################################################################
#............................................................................................................................#
##############################################################################################################################
# FYI, For the limma part, we used a different version of R
# # ml R/4.2.0
# .libPaths("./libraries/R_LIBS_4p2p0/")

#Set/Create Working Directory to Folder
parent_directory <- "/hubmap_single_cell/"
setwd(paste0(parent_directory, "scATAC/projects/"))

# goana in limma
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(future)
library(Matrix)
library(org.Hs.eg.db)
library(limma)
library(circlize)
#library(ggpubr)
library(data.table)
library("AnnotationDbi")
library(BuenColors)

clustering <- readRDS("clustering_genes_traj.rds")

for (i in 1:7){
  cluster_genes <- c(names(which(clustering$cluster == as.character(i))))
  cluster_genes <- transpose(data.frame(strsplit(cluster_genes, ":")))[,2]
  current_clusters <- as.character(i)
  message(length(cluster_genes))

  entrez.ids <- mapIds(org.Hs.eg.db, keys=cluster_genes, 
      column="ENTREZID", keytype="SYMBOL")
  kegg.out <- kegga(unique(entrez.ids), species="Hs")
  kegg.out <- kegg.out[order(kegg.out$P.DE),]
  kegg.useful <- kegg.out[kegg.out$N <= 200,]
  kegg.useful$FDR <- p.adjust(kegg.useful$P.DE, method = "fdr")
  message(head(kegg.useful, 20))

  if (i==1){
    kegg.save <- kegg.useful
  } else {
    kegg.save[,paste0("P.DE", i)] <- kegg.useful[rownames(kegg.save),]$P.DE
  }

  df <- head(kegg.useful, 10)
  df$Pathway <- factor(df$Pathway, levels = df$Pathway)
  df$log10p <- (-log10(df$P.DE))
  df$log10FDR <- (-log10(df$FDR))
  p <- ggplot(df, aes(x = Pathway, y = log10FDR))+
    geom_col(width = 0.7)+ coord_flip()# + theme_ArchR()
  #plotPDF(p, name = paste0("KEGG_pathway_neglog10FDR_cluster", current_clusters,"_", length(cluster_genes), "genes.pdf"), width = 8, height = 5, ArchRProj = proj_duodenum, addDOC = FALSE)
}
top30 <- kegg.save[order(matrixStats::rowMins(as.matrix(kegg.save[,c(4,6:11)])))[1:40],]
rownames(top30) <- top30$Pathway
paletteLength <- 256
# myBreaks <- c(seq(0, 5, length.out=ceiling(paletteLength/2) + 1),
#     seq(5.01, 10, length.out=ceiling(paletteLength/2) + 1))
mat_to_plot <- -log10(as.matrix(top30[,c(4,6:11)]))
mat_to_plot[mat_to_plot>10] <- 10
p <- pheatmap::pheatmap(
    mat = mat_to_plot, 
    color = jdb_palette("solar_extra", type = "continuous"), show_rownames = T, 
    border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE#, breaks = myBreaks
)
ggsave(paste0("kegga","_top40terms", "genes.pdf"), p, width = 8, height = 6)


