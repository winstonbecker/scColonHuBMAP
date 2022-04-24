# scColonHuBMAP  
## Introduction  
Code for snATAC and snRNA analysis in High Resolution Single Cell Maps Reveals Distinct Cell Organization and Function Across Different Regions of the Human Intestine
by Hickey*, Becker*, and Nevins* et al (https://www.biorxiv.org/content/10.1101/2021.11.25.469203v1). This github wiki readme and the corresponding code will walk through the analysis done on the intestine data collected in this paper. 

## Description of dataset and data availability  
This dataset consists of data from 8 regions of the intestine (Duodenum, Proximal-jejunum, Mid-jejunum, Ileum, Ascending Colon, Transverse Colon, Descending Colon, and Sigmoid Colon) from 9 donors (B001, B004, B005, B006, B008, B009, B010, B011, and B012). Matched scATAC and snRNA data were colected from 3 of these donors (B001, B004, and B005) and multiome data (scATAC and scRNA in the same cell) was collected from the remaining 6 donors. The data can be obtained from the HuBMAP data portal and will be made available on dbGaP. Processed data files including expression matricies will be linked to on this page.  

## Quality Control and Filtering 
The first step in analyzing this data was quality control and filtering. For the scATAC data, cells were filtered based on TSS enrichment and number of fragments/droplet. We set a TSS enrichment cutoff of 5 for all samples, and set specific cutoffs for the minimum number of fragments per droplet for each sample due to different sequencing depth in different samples. Different cutoffs were maually determined to isolate a population of cells. Following filtering based on these QC metrics we simulated scATAC doubles using ArchR and removed predicted doublets from the dataset. This anlaysis is done in the script 1_scATAC_initial_QC_and_filtering.R.   

For the snRNA data, we started by first filtering out multiome cells that did not meet the scATAC QC cutoffs. We next filtered cells to have a minimum of *** UMI/droplet.  Finally, we simulated doublets for the non-multiome samples and removed predicted doublets from the datasets. 

## Initial Clustering of scATAC data  

## Initial Clustering of scRNA data

## Subclustering and annotation of scRNA data

## Subclustering and annotation of scATAC data

## Data integration

## Regulatory TF analysis

## Trajectory analysis

## LD score Regression

