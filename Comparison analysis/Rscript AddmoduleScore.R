stringsAsFactors=FALSE
library(Seurat)
library(genefilter)
library(GSEABase)
library(stringr)
library(foreach)
library(GSA)
rm(list=ls())
gc()


#####################################################################################################
############## Reading mRNAExpr into R  ############################################
# load the simulated mRNAExpr
dir_path <- "/gpfs/fs2/scratch/rwang46/Single Cell Pathway Score_scPS/Job 1 Simulation/Practical/Sample size 200/Scenario 1 scimpute/"
setwd(dir_path)

file_list <- list.files(dir_path,recursive = F, full.names = F)
file_exts <- tools::file_ext(file_list)
file_splits <- split(file_list,file_exts)

############################# Reading genes sets and computing scores using SCSE function and dopar  ################################

for (repl in 1:10) {
  load(file = file_splits$RData[repl])
  predBind <- GSEABase::getGmt(file.path(file_splits$gmt[repl]))
  predBind <- geneIds(predBind)
  Participant.integrated.1 <- AddModuleScore(Participant.integrated.1,predBind)
  AddModule_Score <- t(Participant.integrated.1@meta.data[, -which(colnames(Participant.integrated.1@meta.data) %in% c("orig.ident","nCount_RNA","nFeature_RNA","Condition"))])
  rownames(AddModule_Score) <- names(predBind)
  write.csv(AddModule_Score,file = paste("Score/rep",repl,'AddModuleScore.csv',sep = " "))
}  

