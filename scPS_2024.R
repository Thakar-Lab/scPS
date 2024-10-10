## -----------------------------------------------------------------------------
# scPS (single-cell Pathway Score):- 
##   single-cell RNA-seq gene set analysis (scGSA) method

##   Input Data:- 
###   1. A Seurat object (Seurat V5)
###   2. A gene matrix transposed (GMT) file format of gene sets/pathways

##   Output:-
###   1. Single cells with corresponding scores for gene sets/pathways

##     Authors:-  Ruoqiao Wang
##     Email:-    RuoQiao_Wang@URMC.Rochester.edu

## ----How to Apply scPS Function-----------------------------------------------
### Seurat_data = Input Seurat object
### GeneSet = Gene set gmt file
Result  =  scPS(Seurat_data,GeneSet) # calling scPS 

## ----libraries required-------------------------------------------------------
library(Seurat) #Seurat V5
library(GSEABase)
library(genefilter)

## ----scPS functions-----------------------------------------
#### Function 1:- Run Principal Component Analysis of gene set/pathway ####
GS_PCA_Calculation <- function(Seurat_data,GeneSet){
  GS_PCA <- list()
  # Check if scale.data from RNA or SCT assay is present
  if("scale.data" %in% names(Seurat_data@assays$RNA@layers) | "SCT" %in% names(Seurat_data@assays)){
    print("scale.data was applied for PCA")
    # Run Principal Component Analysis for gene set/pathway
    for (currGS in names(GeneSet)){
      GS_PCA[[currGS]] <- RunPCA(Seurat_data, features = geneIds(GeneSet[[currGS]]), npcs = 10, weight.by.var = FALSE, verbose = F) # approx = FALSE
    }
    print("successfully executed PCA")
    # cumulVarExp and CombPC
    varExplThresh <- 0.5 # target for cumulative max var explained
    cumulVarExp <- function(sdevPC){cumsum(sdevPC^2)/sum(sdevPC^2)}
    varExp <- function(sdevPC){(sdevPC^2) / sum(sdevPC^2)}
    maxCombPC <- unlist(lapply(GS_PCA, FUN = function(x){
      idx = min(which(cumulVarExp(x[['pca']]@stdev) > varExplThresh));
      if(identical(idx, integer(0))){return(1)}else{return(max(idx))}
    }))
    # Scale PC score above 0
    aggregPC.scaled <- mat.or.vec(nr = length(names(GS_PCA)),nc = dim(Seurat_data)[2])
    rownames(aggregPC.scaled) <- names(GS_PCA)
    colnames(aggregPC.scaled) <- colnames(Seurat_data)
    PC.X <- list()
    for(currGS in rownames(aggregPC.scaled)){
      PC.X[[currGS]] <- GS_PCA[[currGS]]@reductions$pca@cell.embeddings - min(GS_PCA[[currGS]]@reductions$pca@cell.embeddings)
    }
    # for each gene set after filtering  
    for(currGS in rownames(aggregPC.scaled)){
      # get the var% explained per PC
      currVarExp = varExp(GS_PCA[[currGS]]@reductions$pca@stdev)
      # weight PC scores
      currWeightedPC = sapply(1:maxCombPC[currGS], FUN =
                                function(i){(PC.X[[currGS]][,i] * currVarExp[i])})
      if(ncol(currWeightedPC) > 1){
        aggregPC.scaled[ currGS , rownames(currWeightedPC) ] <-
          sqrt(rowSums(currWeightedPC[ , 1:(maxCombPC[currGS])]))
      }else{
        aggregPC.scaled[ currGS , rownames(currWeightedPC) ] <-
          sqrt(sum(currWeightedPC[ , 1:(maxCombPC[currGS])]))
      }
      rm(currWeightedPC, currVarExp)  
    }} 
  # Error: run ScaleData or SCTransform
  else{
    print("No layer matching pattern 'scale.data' not found. Please run ScaleData or SCTransform and retry")
  }
  return(aggregPC.scaled)
}


#### Function 2:- Calculate the mean expression level of gene set/pathway ####
GS_Experssion_Calculation <- function(Seurat_data,GeneSet,GS_PCs){
  mean_express_list<- list()
  if("SCT" %in% names(Seurat_data@assays)){
    print("normalized data from the SCT was applied")
    for (currGS in rownames(GS_PCs)){
      genes <- geneIds(GeneSet[[currGS]])
      mean_express_list[[currGS]] <- apply(Seurat_data@assays$SCT@data[which(rownames(Seurat_data@assays$SCT@data) %in% genes),],2,mean)
    }
    mean_expression <- mat.or.vec(nr = length(rownames(GS_PCs)),nc = length(colnames(GS_PCs)))
    rownames(mean_expression) <- rownames(GS_PCs)
    colnames(mean_expression) <- colnames(GS_PCs)
    for (currGS in rownames(GS_PCs)){
      mean_expression[currGS,] <- mean_express_list[[currGS]]
    }
    print("successfully calculated the mean expression")
  }else if ("data" %in% names(Seurat_data@assays$RNA@layers)){
    print("normalized data from the RNA was applied")
    for (currGS in rownames(GS_PCs)){
      genes <- geneIds(GeneSet[[currGS]])
      mean_express_list[[currGS]] <- apply(Seurat_data@assays$RNA@layers$data[which(rownames(Seurat_data@assays$RNA@layers$data) %in% genes),],2,mean)
    }
    mean_expression <-mat.or.vec(nr = length(rownames(GS_PCs)),nc = length(colnames(GS_PCs)))
    rownames(mean_expression) <- rownames(GS_PCs)
    colnames(mean_expression) <- colnames(GS_PCs)
    for (currGS in rownames(GS_PCs)){
      mean_expression[currGS,] <- mean_express_list[[currGS]]
    }
    print("successfully calculated the mean expression")
  }else{
    print("No layer matching pattern 'data' not found. Please run NormalizeData or SCTransform and retry")  
  }
  return(mean_expression)
}

#### Function 3:- Calculate scPS score of gene set/pathway ####
scPS <- function(Seurat_data,GeneSet){
  GS_PCs <- GS_PCA_Calculation(Seurat_data,GeneSet)
  GS_MeanExpress <- GS_Experssion_Calculation(Seurat_data,GeneSet,GS_PCs)
  scPS_Scores <- mat.or.vec(nr = length(rownames(GS_PCs)), nc = length(colnames(GS_PCs)))
  rownames(scPS_Scores) <- rownames(GS_PCs)
  colnames(scPS_Scores) <- colnames(GS_PCs)
  for (j in 1:length(colnames(GS_PCs))){
    for (i in 1:length(rownames(GS_PCs))){
      scPS_Scores[i,j] <- GS_PCs[i,j] * GS_MeanExpress[i,j] 
    }
  }
  print("DONE!")
  return(scPS_Scores)
}


# Example here
predBind <- GSEABase::getGmt("KEGG_pathway_2024.gmt")
numTF <- length(unique(names(predBind)))
# FILTER OUT MOTIFS THAT HAVE EITHER TOO MANY OR TOO FEW PREDICTED SITES
predBindCounts <- sort(setNames(unlist(lapply(predBind, function(x){length(geneIds(x))})), names(predBind)))
predBind_maxMin <- list(max = 60, min = 50) # SPECIFY CUTOFFS
predBindCounts_postThresh <- predBindCounts[(predBindCounts >= predBind_maxMin$min)&
                                              (predBindCounts <= predBind_maxMin$max)]
predBind <- predBind[names(predBind) %in% names(predBindCounts_postThresh)]
print(length(predBind))

Result <- scPS(obj,predBind) # calling scPS

