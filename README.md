# scPS (Single-Cell Pathway Score)
<p align="left">
  <img src="scPS.png" width="100">
</p>
scPS is a single-cell RNA-seq gene set analysis (scGSA) method designed to evaluate gene set activity at the single-cell level. It utilizes functions within an R script, eliminating the need for additional installation. The tool integrates seamlessly with Seurat V5 objects for efficient pathway scoring across thousands of cells.

## ✨ Features
	Input:
	•	A Seurat object (Seurat V5 is required)
	•	A gene set file in GMT (Gene Matrix Transposed) format. This file should contain gene set names, descriptions, and the associated genes.
	Output:
	•	A matrix where rows correspond to gene sets/pathways and columns correspond to single cells, with each cell assigned a score for each gene set or pathway.

## 🛠️ How to Use scPS
```R
Seurat_data = Input Seurat object  
GeneSet = Gene set gmt file  
Result = scPS(Seurat_data, GeneSet)  #Calling scPS  
```

## 📚 Required Libraries
```R
library(Seurat)  #Seurat V5
library(GSEABase)
library(genefilter)
```

## 🔗 Data Availability

The raw data for simulation and comparative analysis are available [here](https://drive.google.com/drive/folders/1Gvp4ydnJbHZEDIxLjyt0xrQcbMziwDBF?usp=drive_link)

## 👩‍🔬 Authors

Ruoqiao Wang (Email: RuoQiao_Wang@URMC.Rochester.edu)

## 📄 Citation

:sparkles: **Please cite when using** `scPS` :sparkles:

[Wang, R. H., & Thakar, J. (2024)](https://academic.oup.com/nargab/article/6/3/lqae124/7770961). Comparative analysis of single-cell pathway scoring methods and a novel approach. NAR Genomics and Bioinformatics, 6(3), lqae124


