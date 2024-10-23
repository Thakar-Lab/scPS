# scPS 
### (single-cell Pathway Score) 
#### scPS is a single-cell RNA-seq gene set analysis (scGSA) method using the functions in R script (without installation)

###   Input Data: 
#### 1. A Seurat object (Seurat V5 are required), 
#### 2. A gene matrix transposed (GMT) file format of gene sets/pathways. Each gene set is described by a name, a description, and the genes in the gene set. 

###   Output:
####   1. Single cells with corresponding scores for gene sets/pathways data matrix, where rows correspond to gene sets/pathways and columns correspond to cells.

### How to Apply scPS Function
#### Seurat_data = Input Seurat object
#### GeneSet = Gene set gmt file
Result  =  scPS(Seurat_data,GeneSet) # calling scPS 

### libraries required
library(Seurat) #Seurat V5
library(GSEABase)
library(genefilter)

All the raw data for the data simulation and comparative analysis can be found at the following link: 
https://drive.google.com/drive/folders/1Gvp4ydnJbHZEDIxLjyt0xrQcbMziwDBF?usp=drive_link

###     Authors:Ruoqiao Wang and Juilee Thakar
###     Email:RuoQiao_Wang@URMC.Rochester.edu and juilee_thakar@urmc.rochester.edu

Please cite the following article when using scPS:

Wang, R. H., & Thakar, J. (2024). Comparative analysis of single-cell pathway scoring methods and a novel approach. NAR Genomics and Bioinformatics, 6(3), lqae124.

https://doi.org/10.1093/nargab/lqae124

