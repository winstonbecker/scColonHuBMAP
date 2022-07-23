# scColonHuBMAP  
## Introduction  
Code for snATAC and snRNA analysis in High Resolution Single Cell Maps Reveals Distinct Cell Organization and Function Across Different Regions of the Human Intestine
by [Hickey*, Becker*, and Nevins* et al](https://www.biorxiv.org/content/10.1101/2021.11.25.469203v1). This github wiki readme and the corresponding code will walk through the analysis done on the intestine single-nuclei data collected in this paper.  

## Description of dataset and data availability  
This dataset consists of data from 8 regions of the intestine (Duodenum, Proximal-jejunum, Mid-jejunum, Ileum, Ascending Colon, Transverse Colon, Descending Colon, and Sigmoid Colon) from 9 donors (B001, B004, B005, B006, B008, B009, B010, B011, and B012). Matched snATAC and snRNA data were colected from 3 of these donors (B001, B004, and B005) and multiome data (snATAC and snRNA in the same cell) was collected from the remaining 6 donors. The data can be obtained from the HuBMAP data portal and will be made available on dbGaP. Links for processed data files including expression matricies will included on this page. Samples were initially processed with cellranger and to obtain fragments files for the snATAC data and the raw expression matricies for the snRNA data for downstream analysis.    

## Quality Control and Filtering 
For the scATAC data, cells were filtered based on TSS enrichment and number of fragments/droplet. We set a TSS enrichment cutoff of 5 for all samples, and set specific cutoffs for the minimum number of fragments per droplet for each sample due to different sequencing depth in different samples. Different cutoffs were maually determined to isolate a population of cells. Following filtering based on these QC metrics we simulated scATAC doublets using ArchR and removed predicted doublets from the dataset. This analysis is done in the script scATAC_1a_create_arrows_initial_filtering.R.  

For the snRNA data, we started by first filtering out multiome cells that did not meet the scATAC QC cutoffs. We next filtered all cells to have a minimum of 400 unique genes/droplet. We also removed cells without at least 3 times the number of UMIs as the median number of UMIs in "empty droplets." For example, prior to any filtering the distribution of UMI across all droplets might look something like this:  
<p align="center">
  <img src="https://user-images.githubusercontent.com/15204322/169681077-2865fc5b-1dfd-45bd-99a9-b9074376f07b.png" width="300" height="250">  
</p> 
Where the peak around 250 UMIs/cell would represent counts in empty droplets and the peak around 2K UMIs would like represent counts in cells. In this case we set the minimum number of UMIs to be at least 3x the median number of UMIs in "empty droplets" to minimize ambient RNA contamination, resulting in the following droplets remaining after filtering:  
<p align="center">
  <img src="https://user-images.githubusercontent.com/15204322/169681092-801c2de5-d0a2-44f1-bc9b-3845aa1e8032.png" width="300" height="250"> 
</p> 

Next we simulated potential doublets for the non-multiome samples using [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) and removed predicted doublets from the datasets. This is carried out in 2_scRNA_analysis_per_sample_qc.R.   

## Ambient RNA Correction
Based on the histograms of UMIs/droplet, we suspect that, like most droplet based scRNA datasets, there is likely some ambient RNA contamination in the dataset. Given that removing ambient RNA can be difficult, we approached this in multiple ways, and completed many of our analyses in parallel while utilizing different methods to control for ambient RNA. As outlined above the first step was to simply remove cells without at least 3 times the number of UMIs as the median number of UMIs in "empty droplets." Next we ran multiple ambient RNA detection methods, including [DecontX](https://github.com/campbio/celda) and [SoupX](https://github.com/constantAmateur/SoupX). For SoupX, we campared the results of using the automated contamination estimation for each individual sample and setting the contmination to 20%, 30%, and 40%.

Example (B010-A-405) of soupX correction of MUC2 expression in immune cells (Corrected Raw counts are on the y-axis with raw counts on the x axis. The far left plot is the automated estimate of contamination and the other plots show 20%, 30%, and 40% contaminaiton, increasing from left to right. Points are colored by cell type):  
<img src="https://user-images.githubusercontent.com/15204322/169681217-7b7232cf-917d-4d6e-970f-8172a524af29.png" width="235" height="200">
<img src="https://user-images.githubusercontent.com/15204322/169681234-10c6ee2b-0a91-45fd-930f-f3a982fa8811.png" width="235" height="200">
<img src="https://user-images.githubusercontent.com/15204322/169681246-543f5805-4cb8-4c5d-a708-772365262598.png" width="235" height="200">
<img src="https://user-images.githubusercontent.com/15204322/169681249-b6868b28-341e-4478-8508-7832b1ff14e1.png" width="235" height="200">

Reassuringly, even when setting the contamination fraction to 40%, known marker genes are not typically removed, as shown here for PAX5 at 40% contamination (red represents B cells and Pink represents Plasma cells).  
<p align="center">
  <img src="https://user-images.githubusercontent.com/15204322/169681280-cc2193dd-6a9e-4b31-9848-41eb2b60aadc.png" width="300" height="250">
</p> 

The clustering and annotation of the dataset was largely unaffected by ambient RNA removal, and differential genes were computed both pre and post correction to confirm stability of results. For example, when we compare the results of immune cell annotated using the raw counts and decontx corrected counts, we observe a high degree of concordance between the annotations. 

Based on these initial results, we used the decontx corrected counts for downstream clustering and annotation of the dataset. This is carried out in the 2b and 2c scripts.  


## Initial Clustering of scRNA data
We first aimed to identify cells belonging to the immune, stromal, and epithelial compartments (see 2_rna_analysis_initial_clustering.R). For this first step, we utilized the standard seurat [normalize and scale pipeline](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) and then corrected batch effects between donors using [Harmony](https://github.com/immunogenomics/harmony). We clustered the harmony corrected matricies and assigned each cluster as immune, stromal, or epithelial based on the expression of marker genes. This is carried out in 3a_scRNA_rna_analysis_initial_clustering.R.   

## Subclustering and annotation of scRNA data
### Immune  
Next we subclustered and annotated the immune RNA cells. For the immune cells, we again ran the seurat standard normalize and scale approach and then ran Harmony to account for batch effects. For initial cell annotation, we leveraged labeling the cells with SingleR, labeling the cells with previously published scRNA datasets, comparing markers for each cluster to previously published colon marker sets, examining conanonical marker genes, and labeling the cells with azimuth. In general, all of these approaches gave similar results and we considered all the results wehen assigning the final cluster identites. Notably, following initial clustering there were a few clusters with RNA expression not consistant with them being immune cells, that may represent doublets. We sorted these cells into possible epithelial and possible stromal groups and saved them for downstream analysis with those compartments. This is carried out in 3b_scRNA_subcluster_immune_rna.R.  

### Stromal  
We then subclustered and annotated the stromal RNA cells following a similar strategy, again using thestandard seurat normalize and sclae approach and then running Harmony for batch correction. Similarly, we removed possibly doublet clusters, and for those expressing epithelial markers saved them for downstream clustering with the epithelial datasets. This is carried out in 3c_scRNA_subcluster_stromal_rna.R.   

### Epithelial  
For the epithelial compartment, we started with all cells that were initially designated as epithelial cells as well as cells expressing epithelail markers when we analyzed the stromal and immune compartments. We divided these cells based on whether they were dervied from the duodenum, jejunum, ileum, and colon, and annotated each location seperately. For annotating all of these regions, we took a different approach from the immune and stromal compartments, this time running scTransform on each individual sample and then integrating the individual samples using the seurat functions FindIntegrationAnchors and IntegrateData. This is carried out in 3d-3h.  

#### Enteroendocrine  
After initially annotating enteroendocrine cells in each compartment, we next annotated all enteroendocrine cells together to better identify subtypes of enteroendocrine cells. This is carried out in 3j_scRNA_enteroendocrine.R.  

#### Specialized Secretory  

## Initial Clustering of scATAC data  

Keep only the ATAC cells that we kept the corresponding RNA cells fror the multiome data. Then question is do we cluster everything together or do it seperately. If seperate then we can go ahead and do it, before we have the final list of multiome cells. 

## Subclustering and annotation of scATAC data

## Can you identify more cell types in the colon with multimodaliy data?
Others have reported that in specific cases multimodal date can enable the elucidation of additional cell states. Following our initial clustering and annotation of the data based ont eh single-ome data, we also attempted...

## Data integration
CCA vs dictionary based method. And if dictionary do you map to the same sample as the references using the entire multiome set as a bridge. Probably best to do by compartment based on what we saw for encode, which may be easier anyway. 

## Regulatory TF analysis
We next aimed to identify regulatory TFs for each cell type in each region of the intestine, which is done in the script scATAC_5a_RegulatorsAllRegions.R.  

For the epithelial comparment, we first identified TF regulators individually for each major region of the intestine (duodenum, jejunum, ileum, and colon) as outlined in the [ArchR manual](https://www.archrproject.com/bookdown/identification-of-positive-tf-regulators.html). In this case, we correlated the Gene Expression Matix with the Vierstra Motif Matrix. TF regulators were selected as those with a correlation between these two matricies of 0.5 and a adjusted pvalue of the correlation of less than 0.01 and a motif delta above the 75th percentile. We then included TFs that were regulators in any of the four regions as the final list of TF regulators. We then plot the row z-scores of the Gene expression matrix for the corresponding TFs and the row z-scores of the chrom var deviation z-score for each TF. 

For the immune and stromal compartments, we identified TF regulators once for all sections rather than individually for each section.  


## Trajectory analysis

## LD score Regression

