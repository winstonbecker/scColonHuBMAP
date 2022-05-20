# scColonHuBMAP  
## Introduction  
Code for snATAC and snRNA analysis in High Resolution Single Cell Maps Reveals Distinct Cell Organization and Function Across Different Regions of the Human Intestine
by [Hickey*, Becker*, and Nevins* et al](https://www.biorxiv.org/content/10.1101/2021.11.25.469203v1). This github wiki readme and the corresponding code will walk through the analysis done on the intestine data collected in this paper. 

## Description of dataset and data availability  
This dataset consists of data from 8 regions of the intestine (Duodenum, Proximal-jejunum, Mid-jejunum, Ileum, Ascending Colon, Transverse Colon, Descending Colon, and Sigmoid Colon) from 9 donors (B001, B004, B005, B006, B008, B009, B010, B011, and B012). Matched scATAC and snRNA data were colected from 3 of these donors (B001, B004, and B005) and multiome data (scATAC and scRNA in the same cell) was collected from the remaining 6 donors. The data can be obtained from the HuBMAP data portal and will be made available on dbGaP. Processed data files including expression matricies will be linked to on this page.  

## Quality Control and Filtering 
The first step in analyzing this data was quality control and filtering. For the scATAC data, cells were filtered based on TSS enrichment and number of fragments/droplet. We set a TSS enrichment cutoff of 5 for all samples, and set specific cutoffs for the minimum number of fragments per droplet for each sample due to different sequencing depth in different samples. Different cutoffs were maually determined to isolate a population of cells. Following filtering based on these QC metrics we simulated scATAC doubles using ArchR and removed predicted doublets from the dataset. This analysis is done in the script 1_scATAC_initial_QC_and_filtering.R.  

For the snRNA data, we started by first filtering out multiome cells that did not meet the scATAC QC cutoffs. We next filtered cells to have a minimum of 400 unique genes/droplet. We also removed cells without at least 3 times the number of UMIs as the median number of UMIs in "empty droplets." Finally, we simulated doublets for the non-multiome samples and removed predicted doublets from the datasets.  

Based on the histograms of UMIs/droplet, we suspect there is likely some ambient RNA contamination in the dataset. We approached this in multiple ways. First we removed cells without at least 3 times the number of UMIs as the median number of UMIs in "empty droplets." Next we ran multiple ambient RNA detection methods, including decontX and soupX. For SoupX, we campared the results of using the automated contamination estimation for each individual sample and setting the contmination to 20%, 30%, and 40%. The clustering and annotation of the dataset was largely unaffected by ambient RNA removal, and differential were computed both pre and post correction to confirm stability of results.  
  

## Initial Clustering of scRNA data
We first aimed to identify cells belonging to the immune, stromal, and epithelial compartments. For this first step, we utilized the standard seurat [normalize and scale pipeline](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) and then corrected batch effects between donors using [Harmony](https://github.com/immunogenomics/harmony). We clustered the harmony corrected matricies and assigned each clusted as immune, stromal, or epithelial based on the expression of marker genes.  

## Subclustering and annotation of scRNA data
### Immune  
Next we subclustered and annotated the immune RNA cells.  

### Stromal  

### Epithelial  
#### Duodenum  
#### Jejunum  
#### Ileum  
#### Colon  

## Initial Clustering of scATAC data  

## Subclustering and annotation of scATAC data

## Can you identify more cell types in the colon with multimodaliy data?

## Data integration

## Regulatory TF analysis

## Trajectory analysis

## LD score Regression

