# scColonHuBMAP  
## Introduction  
Code for snATAC and snRNA analysis in High Resolution Single Cell Maps Reveals Distinct Cell Organization and Function Across Different Regions of the Human Intestine
by [Hickey*, Becker*, and Nevins* et al](https://www.biorxiv.org/content/10.1101/2021.11.25.469203v1). This github wiki readme and the corresponding code will walk through the analysis done on the intestine data collected in this paper. 

## Description of dataset and data availability  
This dataset consists of data from 8 regions of the intestine (Duodenum, Proximal-jejunum, Mid-jejunum, Ileum, Ascending Colon, Transverse Colon, Descending Colon, and Sigmoid Colon) from 9 donors (B001, B004, B005, B006, B008, B009, B010, B011, and B012). Matched scATAC and snRNA data were colected from 3 of these donors (B001, B004, and B005) and multiome data (scATAC and scRNA in the same cell) was collected from the remaining 6 donors. The data can be obtained from the HuBMAP data portal and will be made available on dbGaP. Processed data files including expression matricies will be linked to on this page.  

## Quality Control and Filtering 
The first step in analyzing this data was quality control and filtering. For the scATAC data, cells were filtered based on TSS enrichment and number of fragments/droplet. We set a TSS enrichment cutoff of 5 for all samples, and set specific cutoffs for the minimum number of fragments per droplet for each sample due to different sequencing depth in different samples. Different cutoffs were maually determined to isolate a population of cells. Following filtering based on these QC metrics we simulated scATAC doubles using ArchR and removed predicted doublets from the dataset. This analysis is done in the script 1_scATAC_initial_QC_and_filtering.R.  

For the snRNA data, we started by first filtering out multiome cells that did not meet the scATAC QC cutoffs. We next filtered cells to have a minimum of 400 unique genes/droplet. We also removed cells without at least 3 times the number of UMIs as the median number of UMIs in "empty droplets," as shown in the example below: 

Prefilter Example:  
![image](https://user-images.githubusercontent.com/15204322/169681077-2865fc5b-1dfd-45bd-99a9-b9074376f07b.png)  

Postfilter Example:  
![image](https://user-images.githubusercontent.com/15204322/169681092-801c2de5-d0a2-44f1-bc9b-3845aa1e8032.png)  

Next we simulated doublets for the non-multiome samples and removed predicted doublets from the datasets. Based on the histograms of UMIs/droplet, we suspect that like most droplet based scRNA datasets there is likely some ambient RNA contamination in the dataset. Given that removing ambient RNA can be difficult, we approached this in multiple ways, and completed many of our analyses in parallel considered different methods to control for ambient RNA. First we removed cells without at least 3 times the number of UMIs as the median number of UMIs in "empty droplets." Next we ran multiple ambient RNA detection methods, including decontX and soupX. For SoupX, we campared the results of using the automated contamination estimation for each individual sample and setting the contmination to 20%, 30%, and 40%. The clustering and annotation of the dataset was largely unaffected by ambient RNA removal, and differential genes were computed both pre and post correction to confirm stability of results.  

Example (B010-A-405) of soupX correction of MUC2 expression in immune cells (Corrected Raw counts are on the y-axis with raw counts on the x axis. The far left plot is the automated estimate of contamination and the other plots show 20%, 30%, and 40% contaminaiton, increasing from left to right. Points are colored by cell type):  
<img src="https://user-images.githubusercontent.com/15204322/169681217-7b7232cf-917d-4d6e-970f-8172a524af29.png" width="235" height="200">
<img src="https://user-images.githubusercontent.com/15204322/169681234-10c6ee2b-0a91-45fd-930f-f3a982fa8811.png" width="235" height="200">
<img src="https://user-images.githubusercontent.com/15204322/169681246-543f5805-4cb8-4c5d-a708-772365262598.png" width="235" height="200">
<img src="https://user-images.githubusercontent.com/15204322/169681249-b6868b28-341e-4478-8508-7832b1ff14e1.png" width="235" height="200">

Reassuringly, even when setting the contamination fraction to 40%, known marker genes are not typically removed, as shown here for PAX5 at 40% contamination (red represents B cells and Pink represents Plasma cells).
![image](https://user-images.githubusercontent.com/15204322/169681280-cc2193dd-6a9e-4b31-9848-41eb2b60aadc.png)



## Initial Clustering of scRNA data
We first aimed to identify cells belonging to the immune, stromal, and epithelial compartments. For this first step, we utilized the standard seurat [normalize and scale pipeline](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) and then corrected batch effects between donors using [Harmony](https://github.com/immunogenomics/harmony). We clustered the harmony corrected matricies and assigned each clusted as immune, stromal, or epithelial based on the expression of marker genes.  

## Subclustering and annotation of scRNA data
### Immune  
Next we subclustered and annotated the immune RNA cells.  

### Stromal  
We then subclustered and annotated the stromal RNA cells.

### Epithelial  
#### Duodenum  
#### Jejunum  
#### Ileum  
#### Colon  
#### Enteroendocrine  
After initially annotating enteroendocrine cells in each compartment, we next annotated all enteroendocrine cells together to better identify subtypes of enteroendocrine cells. 

## Initial Clustering of scATAC data  

Keep only the ATAC cells that we kept the corresponding RNA cells fror the multiome data. Then question is do we cluster everything together or do it seperately. If seperate then we can go ahead and do it, before we have the final list of multiome cells. 

## Subclustering and annotation of scATAC data

## Can you identify more cell types in the colon with multimodaliy data?
Others have reported that in specific cases multimodal date can enable the elucidation of additional cell states. Following our initial clustering and annotation of the data based ont eh single-ome data, we also attempted...

## Data integration
CCA vs dictionary based method

## Regulatory TF analysis

## Trajectory analysis

## LD score Regression

