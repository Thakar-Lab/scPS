library(Seurat)
library(scImpute)

setwd("/Users/rwang4/Library/CloudStorage/Box-Box/Ruoqiao Thakar Lab/Single-cell Pathway Score (scPS)/Data & results/HIV and B cell/")
load("~/Library/CloudStorage/Box-Box/Ruoqiao Thakar Lab/Single-cell Pathway Score (scPS)/Data & results/HIV and B cell/Naive Normal B cell.RData")
write.csv(mRNA, "Normal_Bcell_count.csv")

scimpute(count_path="/Users/rwang4/Library/CloudStorage/Box-Box/Ruoqiao Thakar Lab/Single-cell Pathway Score (scPS)/Data & results/HIV and B cell/Normal_Bcell_count.csv",
         infile = "csv",
         outfile = "csv",
         out_dir = "./",
         Kcluster = 1,
         drop_thre = 0.5,
         ncores = 4)