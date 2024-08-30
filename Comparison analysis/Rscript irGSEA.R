## ---- Load package-------------------------------------------------
library(gage)
library(irGSEA)
library(Seurat)
library(genefilter)
library(GSEABase)
rm(list=ls()) 
gc()

## ----Assign gene signature for simulation--------------------------------------------------
# load the simulated data
dir_path <- "/gpfs/fs2/scratch/rwang46/Single Cell Pathway Score_scPS/Job 1 Simulation/Practical/Sample size 200/Scenario 2 scimpute/"
setwd(dir_path)
dir.create("Score")

file_list <- list.files(dir_path,recursive = F, full.names = F)
file_exts <- tools::file_ext(file_list)
file_splits <- split(file_list,file_exts)

for (repl in 1:10){
  load(file = file_splits$RData[repl])
  predBind <- GSEABase::getGmt(file.path(file_splits$gmt[repl]))
  predBind <- geneIds(predBind)
  Participant.integrated.1 <- irGSEA.score(object = Participant.integrated.1, assay = "RNA", slot = "data",
                                           seeds = 42,
                                           custom = T, geneset = predBind,
                                           ncores = 8,
                                           method = c("AUCell","JASMINE","ssgsea", "UCell","scSE"),
                                           kcdf = 'Gaussian',
                                           ucell.MaxRank = NULL,
                                           aucell.MaxRank = NULL,
                                           JASMINE.method = "oddsratio")


  Seurat::Assays(Participant.integrated.1)
  write.csv(Participant.integrated.1@assays$scSE@counts,file = paste("Score/rep",repl,'scSE.csv',sep = " "))
  write.csv(Participant.integrated.1@assays$UCell@counts,file = paste("Score/rep",repl,'UCell.csv',sep = " "))
  write.csv(Participant.integrated.1@assays$ssgsea@counts,file = paste("Score/rep",repl,'ssGSEA.csv',sep = " "))
  write.csv(Participant.integrated.1@assays$JASMINE@counts,file = paste("Score/rep",repl,'JASMINE.csv',sep = " "))
  write.csv(Participant.integrated.1@assays$AUCell@counts,file = paste("Score/rep",repl,'AUCell.csv',sep = " "))
}



