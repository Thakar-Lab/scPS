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
  print("ScaleData is applied using all features")
  # Seurat_data <- ScaleData(Seurat_data, features = rownames(Seurat_data))
  # Check if scale.data from RNA or SCT assay is present
  if("scale.data" %in% names(Seurat_data@assays$RNA@layers) | "SCT" %in% names(Seurat_data@assays)){
    print("scale.data was applied for PCA")
    if(length(GeneSet)>1){
      GSdiscription <- names(GeneSet)
    }else{
      GSdiscription <- setName(GeneSet)
    }
    # Run Principal Component Analysis for gene set/pathway
    for (currGS in GSdiscription){
      if(length(GeneSet)>1){
        genes <- geneIds(GeneSet[[currGS]])
      }else{
        genes <- geneIds(GeneSet)
      }
        print(paste("executed PCA for",currGS))
        GS_PCA[[currGS]] <- RunPCA(Seurat_data, features = genes, npcs = 10, weight.by.var = FALSE, verbose = F) # approx = FALSE
        print(paste("successfully executed PCA for",currGS))
      }
    # cumulVarExp and CombPC
    varExplThresh <- 0.5 # target for cumulative max var explained
    cumulVarExp <- function(sdevPC){cumsum(sdevPC^2)/sum(sdevPC^2)}
    varExp <- function(sdevPC){(sdevPC^2) / sum(sdevPC^2)}
    maxCombPC <- unlist(lapply(GS_PCA, FUN = function(x){
      idx = min(which(cumulVarExp(x[['pca']]@stdev) > varExplThresh));
      if(identical(idx, integer(0))){return(1)}else{return(max(idx))}
    }))
    # Scale PC score above 0
    aggregPC.scaled <- mat.or.vec(nr = length(GSdiscription),nc = dim(Seurat_data)[2])
    rownames(aggregPC.scaled) <- GSdiscription
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
          sqrt(currWeightedPC[ , 1:(maxCombPC[currGS])])
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
GS_Experssion_Calculation <- function(Seurat_data,GeneSet){
  mean_express_list<- list()
  if("SCT" %in% names(Seurat_data@assays)){
    print("normalized data from the SCT was applied")
    if(length(GeneSet)>1){
      GSdiscription <- names(GeneSet)
    }else{
      GSdiscription <- setName(GeneSet)
    }
    for (currGS in GSdiscription){
      if(length(GeneSet)>1){
        genes <- geneIds(GeneSet[[currGS]])
      }else{
        genes <- geneIds(GeneSet)
      }
      mean_express_list[[currGS]] <- apply(Seurat_data@assays$SCT@data[which(rownames(Seurat_data@assays$SCT@data) %in% genes),],2,mean)
    }
    mean_expression <- mat.or.vec(nr = length(GSdiscription),nc = length(colnames(Seurat_data)))
    rownames(mean_expression) <- GSdiscription
    colnames(mean_expression) <- colnames(Seurat_data)
    for (currGS in GSdiscription){
      mean_expression[currGS,] <- mean_express_list[[currGS]]
    }
    print("successfully calculated the mean expression")
  }else if ("data" %in% names(Seurat_data@assays$RNA@layers)){
    print("normalized data from the RNA was applied")
    if(length(GeneSet)>1){
      GSdiscription <- names(GeneSet)
    }else{
      GSdiscription <- setName(GeneSet)
    }
    for (currGS in GSdiscription){
      if(length(GeneSet)>1){
        genes <- geneIds(GeneSet[[currGS]])
      }else{
        genes <- geneIds(GeneSet)
      }
      mean_express_list[[currGS]] <- apply(Seurat_data@assays$RNA@layers$data[which(rownames(Seurat_data@assays$RNA@layers$data) %in% genes),],2,mean)
    }
    mean_expression <-mat.or.vec(nr = length(GSdiscription),nc = length(colnames(Seurat_data)))
    rownames(mean_expression) <- GSdiscription
    colnames(mean_expression) <- colnames(Seurat_data)
    for (currGS in GSdiscription){
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
  GS_MeanExpress <- GS_Experssion_Calculation(Seurat_data,GeneSet)
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

# load KEGG pathway
predBind <- GSEABase::getGmt("./GSA/KEGG_pathway_2024.gmt")
numTF <- length(unique(names(predBind)))
# FILTER OUT MOTIFS THAT HAVE EITHER TOO MANY OR TOO FEW PREDICTED SITES
predBindCounts <- sort(setNames(unlist(lapply(predBind, function(x){length(geneIds(x))})), names(predBind)))
predBind_maxMin <- list(max = 200, min = 190) # SPECIFY CUTOFFS
predBindCounts_postThresh <- predBindCounts[(predBindCounts >= predBind_maxMin$min)&
                                              (predBindCounts <= predBind_maxMin$max)]
predBind <- predBind[names(predBind) %in% names(predBindCounts_postThresh)]
print(length(predBind))

##### scPS ############################################
Result <- scPS(obj,predBind) # calling scPS

# output as csv file
write.csv(predBind_score,file ='predBind_score.csv')
# save it as a new assay in seuratobj
obj[["scPS"]] <- CreateAssayObject(data = predBind_score)

##### note ############################################
# if you get the NA value, please check your rownames(Seurat_data[["RNA"]]@layers$data) or rownames(Seurat_data[["SCT"]]@layers$data) != NULL
# or assign the data with the correct gene name
# rownames(Seurat_data[["RNA"]]@layers$data) = rownames(obj[["RNA"]]) 
DefaultAssay(obj) <- "RNA"
rownames(obj@assays$RNA@layers$data) <- rownames(obj)
rownames(obj@assays$SCT@layers$data) <- rownames(obj)


# if you have trouble with some weired gene name in gene set just remove it from the list
# Optional 2: remove from the list
selected_pathways <- c("hsa04062_Chemokine_signaling_pathway", "hsa04613_Neutrophil_extracellular_trap_formation")
predBind <- predBind[!names(predBind) %in% selected_pathways]
# Final count
print(length(predBind))

# Optional 2: filter specific pathways only
selected_pathways <- c("hsa04062_Chemokine_signaling_pathway", "hsa04613_Neutrophil_extracellular_trap_formation")
predBind <- predBind[names(predBind) %in% selected_pathways]
# Final count
print(length(predBind))



