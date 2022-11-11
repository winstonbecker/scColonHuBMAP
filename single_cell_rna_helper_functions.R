# WRB 2022
##############################################################################################################################
# Function Definitions
seurat_standard_normalize_and_scale <- function(colon, cluster = TRUE, cluster_resolution = 1.0, n_dims = 20){
	# Function to run standard seurat pipeline
	# colon: seurat object to run through standard pipeline

	colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
	colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(colon)
	colon <- ScaleData(colon, features = all.genes)
	colon <- RunPCA(colon, features = VariableFeatures(object = colon))
	if (cluster){
		colon <- FindNeighbors(colon, dims = 1:n_dims)
		colon <- FindClusters(colon, resolution = cluster_resolution)
	}
	colon <- RunUMAP(colon, dims = 1:n_dims)
	return(colon)
}

make_seurat_object_qc_only <- function(colon.data, project_name){
	# Function for seurat object creation and basic qc
	# data_directory: contains the expression matrix
	# project_name: name that will be used for project in the seurat object and that will be used when saving plots

	currentSample <- CreateSeuratObject(counts = colon.data, project = project_name, min.cells = 3, min.features = 1)
	currentSample[["percent.mt"]] <- PercentageFeatureSet(currentSample, pattern = "^MT-")

	# keep top 100K barcodes
	if (dim(currentSample)[2]>100000){
		currentSample <- subset(currentSample, subset = nCount_RNA > currentSample$nCount_RNA[rev(order(currentSample$nCount_RNA))][100000])
	}
	
	# Plot some QC plots priot to filtering: QC violin plots, feature scatter, feature histogram, count histogram, barcode plots
	pdf(paste0("./", project_name, "qc_plots", "_prefiltered_no_points.pdf"))
	print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
	dev.off()
	pdf(paste0("./", project_name, "qc_plots", "_feature_scatter.pdf"))
	print(FeatureScatter(currentSample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
	dev.off()
	pdf(paste0("./", project_name, "qc_plots", "_initial_nCount_hist.pdf"))
	hist(log10(currentSample$nCount_RNA+1),main='counts per cell',col='#272E6A')
	dev.off()
	pdf(paste0("./", project_name, "qc_plots", "_initial_nFeature_hist.pdf"))
	hist(log10(currentSample$nFeature_RNA+1),main='genes per cell',col='#272E6A')
	dev.off()
	# currentSample <- CalculateBarcodeInflections(currentSample)
	# pdf(paste0("./", project_name, "BarcodeInflectionsPlot", ".pdf"))
	# BarcodeInflectionsPlot(currentSample)
	# dev.off()
	# currentSample <- CalculateBarcodeInflections(currentSample, threshold.high = 10000)
	# pdf(paste0("./", project_name, "BarcodeInflectionsPlot_10K_threshold", ".pdf"))
	# BarcodeInflectionsPlot(currentSample)
	# dev.off()

	# Now filter everything to greater than 3x the median of the counts in the top 100K droplets, great than 400 unique genes, mt% less than 5, max features of 10K and max counts of 20K
	initial_med <- median(currentSample$nCount_RNA)
	currentSample <- subset(currentSample, subset = nCount_RNA > initial_med*3 & nFeature_RNA > 400 & nFeature_RNA < 10000 & percent.mt < 5 & nCount_RNA < 20000)

	message(paste0("Sample filtered to at least ", initial_med*3, " counts and at least 400 features."))

	# Plot histograms of the results
	pdf(paste0("./", project_name, "qc_plots", "_3x_med_counts_nCount_hist.pdf"))
	hist(log10(currentSample$nCount_RNA+1),main='counts per cell',col='#272E6A')
	dev.off()
	pdf(paste0("./", project_name, "qc_plots", "_3x_med_counts_nFeature_hist.pdf"))
	hist(log10(currentSample$nFeature_RNA+1),main='genes per cell',col='#272E6A')
	dev.off()

	# Return the filtered seurat object
	return(currentSample)
}

make_seurat_object_and_doublet_removal <- function(data_directory, project_name, runDoubletFinder = FALSE, runSoupX = FALSE, cell_filter = NULL, passed_immune_labels = NULL, gene_to_check = c("PAX5", "TMPRSS15", "MUC2", "MUC6", "ACTA2", "CD8A", "CLCA4", "SLT1", "ALDO3", "ROR1")){
	# Function for seurat object creation and basic qc and filtering and running of doublet finder
	# data_directory: contains the expression matrix
	# project_name: name that will be used for project in the seurat object and that will be used when saving plots

	# load data
	colon.data <- Read10X(data.dir = data_directory)
	if (length(colon.data) == 2){
		colon.data <- colon.data[["Gene Expression"]]
	}

	# seurat doesn't like underscores in feature names
	library(stringr)
	rownames(colon.data) <- str_replace_all(rownames(colon.data), "_", "-")
	
	# Make seurat object and filter
	currentSample <- make_seurat_object_qc_only(colon.data, project_name)
	assays_to_keep <- "RNA"

	if (!is.null(cell_filter)){
		message(paste0(dim(currentSample)[2], " cells prior to filtering by list"))
		currentSample <- currentSample[,colnames(currentSample) %in% cell_filter]
		message(paste0(dim(currentSample)[2], " cells remaining following filtering by list"))
	}
	if (runDoubletFinder){
		# Standard normalization and UMAP in preperation for running doublet finder
		currentSample <- seurat_standard_normalize_and_scale(currentSample, cluster = FALSE)

		# Run doublet finder
		set.seed(1)
		sweep.res.list <- paramSweep_v3(currentSample, PCs = 1:20, sct = FALSE)
		sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
		bcmvn <- find.pK(sweep.stats)
		pK <- as.numeric(as.character(bcmvn$pK))[which.max(bcmvn$BCmetric)]
		nExp_poi <- round(0.076*length(currentSample@meta.data$orig.ident)*length(currentSample@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
		currentSample <- doubletFinder_v3(currentSample, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
		print(head(currentSample@meta.data))
		
		# Rename columns to simple names that can be combined more easily
		currentSample$doublet.class <- currentSample[[paste0("DF.classifications_0.25_", pK, "_",nExp_poi)]]
		currentSample[[paste0("DF.classifications_0.25_", pK, "_",nExp_poi)]] <- NULL
		pann <- grep(pattern="^pANN", x=names(currentSample@meta.data), value=TRUE)
		currentSample$pANN <- currentSample[[pann]]
		currentSample[[pann]] <- NULL

		# Plot pre and post doublet finder results - inspect to make sure results look reasonable
		pdf(paste0("./", project_name, "_UMAP_pre_double_removal", ".pdf"))
		print(DimPlot(currentSample, reduction = "umap", group.by = "doublet.class", cols = c("#D51F26", "#272E6A")))
		dev.off()
		currentSample <- subset(currentSample, subset = doublet.class != "Doublet")
		pdf(paste0("./", project_name, "_UMAP_post_double_removal", ".pdf"))
		print(DimPlot(currentSample, reduction = "umap", cols = c("#D51F26")))
		dev.off()
	}

	saveRDS(currentSample@assays$RNA@counts, paste0("filtered_counts_matrix", project_name, ".rds"))

	# if (runDecontX){
	# 	filtered_counts_matrix <- currentSample@assays$RNA@counts
	# 	raw_counts_matrix <- colon.data[rownames(toc),]
	# 	decontXresults <- decontX(filtered_counts_matrix, background = raw_counts_matrix)
	# }

	if (runSoupX){
		# Run SoupX
		library(SoupX)
		currentSample <- seurat_standard_normalize_and_scale(currentSample)#, cluster_resolution = 2)
		pdf(paste0("./", project_name, "_UMAP_pre_soupX", ".pdf"))
		print(DimPlot(currentSample, reduction = "umap"))
		dev.off()
		toc <- currentSample@assays$RNA@counts
		tod <- colon.data[rownames(toc),]
		sc <- SoupChannel(tod, toc)
		sc <- setClusters(sc, setNames(currentSample$seurat_clusters, names(currentSample$seurat_clusters)))

		# Soup X
		#sc <- autoEstCont(sc, forceAccept = TRUE)#,  tfidfMin = 0.5)
		an.error.occured <- FALSE
		pdf(paste0("./", project_name, "_rho_est", ".pdf"))
		tryCatch( { sc = autoEstCont(sc, priorRho = 0.15)}
	          , error = function(e) {an.error.occured <<- TRUE})
		dev.off()
		if (an.error.occured){
			nCounts <- colSums(colon.data)
			amb_med <- median(nCounts[rev(order(nCounts))][10000:100000])
			cont_med <- median(currentSample$nCount_RNA)
			sc <- setContaminationFraction(sc, min(amb_med/cont_med, 0.3))
			message(paste0("Contamination Fraction Forced to ", min(amb_med/cont_med, 0.3)))
		}
		out <- adjustCounts(sc)

		# create a new assay to store correct
		soupXcounts <- CreateAssayObject(counts = out)

		# add this assay to the previously created Seurat object
		currentSample[["soupXcounts"]] <- soupXcounts

		# now do with predefined contamination fractions that are higher
		toc <- currentSample@assays$RNA@counts # filtered cells
		tod <- colon.data[rownames(toc),] # all cells, keep only the genes in the filtered set
		sc <- SoupChannel(tod, toc)
		sc <- setClusters(sc, setNames(currentSample$seurat_clusters, names(currentSample$seurat_clusters)))
		sc <- setContaminationFraction(sc, 0.2)
		out <- adjustCounts(sc)
		soupXcounts <- CreateAssayObject(counts = out)
		currentSample[["soupXcounts0p2"]] <- soupXcounts

		toc <- currentSample@assays$RNA@counts # filtered cells
		tod <- colon.data[rownames(toc),] # all cells, keep only the genes in the filtered set
		sc <- SoupChannel(tod, toc)
		sc <- setClusters(sc, setNames(currentSample$seurat_clusters, names(currentSample$seurat_clusters)))
		sc <- setContaminationFraction(sc, 0.3)
		out <- adjustCounts(sc)
		soupXcounts <- CreateAssayObject(counts = out)
		currentSample[["soupXcounts0p3"]] <- soupXcounts

		toc <- currentSample@assays$RNA@counts # filtered cells
		tod <- colon.data[rownames(toc),] # all cells, keep only the genes in the filtered set
		sc <- SoupChannel(tod, toc)
		sc <- setClusters(sc, setNames(currentSample$seurat_clusters, names(currentSample$seurat_clusters)))
		sc <- setContaminationFraction(sc, 0.4)
		out <- adjustCounts(sc)
		soupXcounts <- CreateAssayObject(counts = out)
		currentSample[["soupXcounts0p4"]] <- soupXcounts

		assays_to_keep <- c("RNA", "soupXcounts", "soupXcounts0p2", "soupXcounts0p3", "soupXcounts0p4")

		if (!is.null(passed_immune_labels)){
			immune_cells <- passed_immune_labels[grepl(paste0(project_name, "_"), rownames(passed_immune_labels)),, drop = FALSE]
			keep <- str_split_fixed(rownames(immune_cells), paste0(project_name, "_"),2)[,2]

			if (length(keep)>0){
				for (gene in gene_to_check){
					if (gene %in% rownames(currentSample)){
						# note adding a little scatter to the x axis so its easier to see how many points are there
						df <- data.frame(RNA = (GetAssayData(object = currentSample, slot = "counts", assay = "RNA")[gene,keep]),
							soupx = GetAssayData(object = currentSample, slot = "counts", assay = "soupXcounts")[gene,keep],
							soupx0p2 = GetAssayData(object = currentSample, slot = "counts", assay = "soupXcounts0p2")[gene,keep],
							soupx0p3 = GetAssayData(object = currentSample, slot = "counts", assay = "soupXcounts0p3")[gene,keep],
							soupx0p4 = GetAssayData(object = currentSample, slot = "counts", assay = "soupXcounts0p4")[gene,keep],
							celltype = immune_cells$CellType)

						df$RNA <- df$RNA + rnorm(length(keep), mean = 0, sd = 0.05)

						p <- ggplot(df, aes(x=RNA, y=soupx0p2, color = celltype)) + geom_point() + xlim(-1,max(c(df$RNA, 1))) + ylim(-1,max(c(df$RNA, 1))) + theme_ArchR() + geom_abline(slope=1, intercept=0) + scale_color_manual(values = paste0(ArchRPalettes$stallion))
						ggsave(paste0("./", project_name, "_", gene, "_RNA_soupx0p2_scatter_", length(keep) ,"cells.pdf"), plot = p, width = 4, height = 4, useDingbats=FALSE)

						p <- ggplot(df, aes(x=RNA, y=soupx0p3, color = celltype)) + geom_point() + xlim(-1,max(c(df$RNA, 1))) + ylim(-1,max(c(df$RNA, 1))) + theme_ArchR() + geom_abline(slope=1, intercept=0) + scale_color_manual(values = paste0(ArchRPalettes$stallion))
						ggsave(paste0("./", project_name, "_", gene, "_RNA_soupx0p3_scatter_", length(keep) ,"cells.pdf"), plot = p, width = 4, height = 4, useDingbats=FALSE)

						p <- ggplot(df, aes(x=RNA, y=soupx0p4, color = celltype)) + geom_point() + xlim(-1,max(c(df$RNA, 1))) + ylim(-1,max(c(df$RNA, 1))) + theme_ArchR() + geom_abline(slope=1, intercept=0) + scale_color_manual(values = paste0(ArchRPalettes$stallion))
						ggsave(paste0("./", project_name, "_", gene, "_RNA_soupx0p4_scatter_", length(keep) ,"cells.pdf"), plot = p, width = 4, height = 4, useDingbats=FALSE)

						p <- ggplot(df, aes(x=RNA, y=soupx, color = celltype)) + geom_point() + xlim(-1,max(c(df$RNA, 1))) + ylim(-1,max(c(df$RNA, 1))) + theme_ArchR() + geom_abline(slope=1, intercept=0) + scale_color_manual(values = paste0(ArchRPalettes$stallion))
						ggsave(paste0("./", project_name, "_", gene, "_RNA_soupx_scatter_", length(keep) ,"cells.pdf"), plot = p, width = 4, height = 4, useDingbats=FALSE)
					}
				}
			}
		}
	}

	# Remove extra stuff and return filtered Seurat object
	currentSample <- DietSeurat(currentSample, counts=TRUE, data=TRUE, scale.data=FALSE, assays=assays_to_keep)
	return(currentSample)
}

seurat_feature_plot <- function(colon, sample_name, reduction, cell_type, markers){
	# Function to organize seurat feature plots and plot the right size figure when plotting a bunch of feature plots at once
	p1 <- FeaturePlot(colon, features = markers, reduction = reduction, sort.cell = TRUE, combine = FALSE, pt.size = 2)
	fix.sc <- scale_colour_gradientn(colours = paletteContinuous(set = "blueYellow", n = 256, reverse = FALSE))
	if (length(p1)==1){
		width <- 4
		height <- 4
	} else if (length(p1)==2){
		width <- 8
		height <- 4
	} else if (length(p1)<5){
		width <- 8
		height <- 8
	} else if (length(p1)<7){
		width <- 12
		height <- 8
	} else if (length(p1)<10){
		width <- 12
		height <- 12
	} else if (length(p1)<13){
		width <- 16
		height <- 12
	} else if (length(p1)<17){
		width <- 16
		height <- 16
	}
	pdf(paste0("./", reduction, "_feature_plot_", sample_name, "_", cell_type ,".pdf"), width = width, height = height)
	print(CombinePlots(lapply(p1, function (x) AugmentPlot(x + fix.sc))))
	dev.off()
}

nice_qc_violin_plots <- function(colon, sample_name){
	# Function to make some qc plots
	# colon: seurat object
	# sample_name: name for saving figures
	pdf(paste0("./n_genes_violin_", sample_name, ".pdf"), width = 6, onefile=F)
	print(VlnPlot(colon, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",length(unique(colon@meta.data$orig.ident)))))+
	geom_boxplot(outlier.shape = NA, alpha = 0.6)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
	scale_x_discrete(labels=paste0(data.frame(table(colon@meta.data$orig.ident))$Var1, "\n n = ",  data.frame(table(colon@meta.data$orig.ident))$Freq)))
	dev.off()
	pdf(paste0("./n_counts_violin_", sample_name, ".pdf"), width = 6, onefile=F)
	print(VlnPlot(colon, features = c("nCount_RNA"), group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",length(unique(colon@meta.data$orig.ident)))))+
	geom_boxplot(outlier.shape = NA, alpha = 0.6)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
	scale_x_discrete(labels=paste0(data.frame(table(colon@meta.data$orig.ident))$Var1, "\n n = ",  data.frame(table(colon@meta.data$orig.ident))$Freq)))
	dev.off()
	pdf(paste0("./pMT_violin_", sample_name, ".pdf"), width = 6, onefile=F)
	print(VlnPlot(colon, features = c("percent.mt"), group.by = "orig.ident", pt.size = 0, cols = c(rep("#D51F26",length(unique(colon@meta.data$orig.ident)))))+
	geom_boxplot(outlier.shape = NA, alpha = 0.6)+theme_ArchR()+theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1))+
	scale_x_discrete(labels=paste0(data.frame(table(colon@meta.data$orig.ident))$Var1, "\n n = ",  data.frame(table(colon@meta.data$orig.ident))$Freq)))
	dev.off()
}

sc_transform_umap_cluster <- function (colon){
	colon[["percent.ribo"]] <- PercentageFeatureSet(colon, pattern = "^RP[SL]")
	colon[["percent.mt"]] <- PercentageFeatureSet(colon, pattern = "^MT-")
	colon <- SCTransform(colon, method = "glmGamPoi", vars.to.regress = c("percent.mt", "percent.ribo"), verbose = TRUE)
	colon <- RunPCA(colon, verbose = TRUE)
	colon <- RunUMAP(colon, dims = 1:30, verbose = TRUE)
	colon <- FindNeighbors(colon, dims = 1:30, verbose = TRUE)
	colon <- FindClusters(colon, verbose = TRUE)
	return(colon)
}

run_harmony <- function (colon, gBy = "orig.ident", resolution = 1.0){
	colon <- RunHarmony(colon, gBy) # will use PCA
	colon <- RunUMAP(colon, dims = 1:20, reduction = "harmony", reduction.name = "umapharmony")
	colon <- FindNeighbors(colon, reduction = "harmony", dims = 1:20)
	colon <- FindClusters(colon, resolution = resolution)
	return(colon)
}

plotUMAP <- function(colon, reduction = "umap", save_name = NULL){
	pdf(paste0("./", reduction, "_clustering_", save_name, ".pdf"))
	plot = DimPlot(colon, reduction = reduction,  group.by = "seurat_clusters", cols = paletteDiscrete(values = unique(colon@meta.data$seurat_clusters), set = "stallion", reverse = FALSE))
	plot = (LabelClusters(plot = plot, id = "seurat_clusters"))
	print(plot)
	dev.off()

	pdf(paste0("./", reduction, "_samples_", save_name, ".pdf"), width = 12)
	print(DimPlot(colon, reduction = reduction, group.by = "orig.ident", cols = paletteDiscrete(values = unique(colon@meta.data$orig.ident), set = "stallion", reverse = FALSE)))
	dev.off()
}

run_singleR_and_plot <- function(colon, reduction = "umapharmony", types_to_use){
	ref <- HumanPrimaryCellAtlasData()
	ref <- ref[,(colData(ref)$label.main %in% types_to_use)]
	singler.pred <- SingleR(test = as.SingleCellExperiment(colon), ref = ref, labels = ref$label.fine)
	colon <- AddMetaData(colon, metadata = singler.pred$labels, col.name = "SingleR.labels")

	pdf(paste0("./", reduction, "_cell_type_singleR_unlabeled.pdf"), width = 20)
	plot = DimPlot(colon, reduction = reduction, group.by = "SingleR.labels", pt.size = 1)	
	print(plot)
	dev.off()

	colon_test <- colon[,colon@meta.data$SingleR.labels %in% names(which(table(singler.pred$labels)>100))]
	pdf(paste0("./", reduction, "_cell_type_singleR_unlabeled20.pdf"), width = 12)
	plot = DimPlot(colon_test, reduction = reduction, group.by = "SingleR.labels", pt.size = 1)
	print(plot)
	dev.off()

	return(colon)
}

find_anchors_label_transfer <- function(colon, seRNA_path, label, reduction = "umapharmony", pred_name = "CellType"){
	seRNA <- readRDS(seRNA_path)
	seRNA[["nCount_RNA"]] = colSums(x = seRNA, slot = "counts")  # nCount_RNA
	seRNA[["nFeature_RNA"]] = colSums(x = GetAssayData(object = seRNA, slot = "counts") > 0)
	seRNA[["percent.mt"]] <- PercentageFeatureSet(seRNA, pattern = "^MT-")
	seRNA <- SCTransform(seRNA, method = "glmGamPoi", vars.to.regress = c("percent.mt"), verbose = FALSE)
	seRNA <- RunPCA(seRNA, verbose = FALSE)
	seRNA <- FindVariableFeatures(object = seRNA)
	colon.anchors <- FindTransferAnchors(reference = seRNA, query = colon, dims = 1:30, normalization.method = 'SCT')
	if (pred_name == "CellType"){
		predictions <- TransferData(anchorset = colon.anchors, refdata = seRNA$CellType, dims = 1:30)
	} else if (pred_name == "Cluster") {
		predictions <- TransferData(anchorset = colon.anchors, refdata = seRNA$Cluster, dims = 1:30)
	}
	colnames(predictions) <- paste0(colnames(predictions), ".", label)
	colon <- AddMetaData(colon, metadata = predictions[,c(paste0("predicted.id.", label)), drop = FALSE])

	pdf(paste0("./UMAP_predicted_id_label_transfer-", label, ".pdf"), width = 12)
	print(DimPlot(colon, reduction = reduction, group.by = paste0("predicted.id.", label), cols = paletteDiscrete(values = unique(colon@meta.data[,paste0("predicted.id.", label)]), set = "stallion", reverse = FALSE)))
	dev.off()

	return(colon)
}


