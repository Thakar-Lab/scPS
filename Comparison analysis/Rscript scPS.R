#### scRNAseq datasets and gene signature simulation ####
## ---- Load package------------------------------------------------------------
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

for (repl in 1:10) {
    load(file = file_splits$RData[repl])
    predBind <- GSEABase::getGmt(file.path(file_splits$gmt[repl]))
    numTF <- length(unique(names(predBind)))
    # FILTER OUT MOTIFS THAT HAVE EITHER TOO MANY OR TOO FEW PREDICTED SITES
    predBindCounts <- sort(setNames(unlist(lapply(predBind, function(x){length(geneIds(x))})), names(predBind)))
    predBind_maxMin <- list(max = 10000, min = 10) # SPECIFY CUTOFFS
    predBindCounts_postThresh <- predBindCounts[(predBindCounts >= predBind_maxMin$min)&
                                                  (predBindCounts <= predBind_maxMin$max)]
    predBind <- predBind[names(predBind) %in% names(predBindCounts_postThresh)]
    print(length(predBind))
    
    mRNAExpr <- Participant.integrated.1@assays$RNA@data
    
    tfTgtPca<- list()
    # RunPCA (with principal component = 10)
    tryCatch(
      {for (currTF in names(predBind)){
        tfTgtPca[[currTF]] <- RunPCA(Participant.integrated.1, features = geneIds(predBind[[currTF]]), npcs = 10, weight.by.var = FALSE, verbose = FALSE) # approx = FALSE
      }
        print("Successfully executed PCA")
        # cumulVarExp and CombPC
        varExplThresh <- 0.5 # target for cumulative max var explained
        cumulVarExp <- function(sdevPC){cumsum(sdevPC^2)/sum(sdevPC^2)}
        varExp <- function(sdevPC){(sdevPC^2) / sum(sdevPC^2)}
        maxCombPC <- unlist(lapply(tfTgtPca, FUN = function(x){
          idx = min(which(cumulVarExp(x[['pca']]@stdev) > varExplThresh));
          if(identical(idx, integer(0))){return(1)}else{return(max(idx))}
        }))
        table(maxCombPC)
        # Scale PC score above 0
        aggregPC.unscaled <- mat.or.vec(nr = length(names(tfTgtPca)),nc = ncol(mRNAExpr))
        rownames(aggregPC.unscaled) <- names(tfTgtPca)
        colnames(aggregPC.unscaled) <- colnames(mRNAExpr)
        PC.X <- list()
        for(currTF in rownames(aggregPC.unscaled)){
          PC.X[[currTF]] <- tfTgtPca[[currTF]]@reductions$pca@cell.embeddings - min(tfTgtPca[[currTF]]@reductions$pca@cell.embeddings)
        }
        # for each gene set after filtering  
        for(currTF in rownames(aggregPC.unscaled)){
          # get the var% explained per PC
          currVarExp = varExp(tfTgtPca[[currTF]]@reductions$pca@stdev)
          # weight PC scores
          currWeightedPC = sapply(1:maxCombPC[currTF], FUN =
                                    function(i){(PC.X[[currTF]][,i] * currVarExp[i])})
          if(ncol(currWeightedPC) > 1){
            aggregPC.unscaled[ currTF , rownames(currWeightedPC) ] <-
              sqrt(rowSums(currWeightedPC[ , 1:(maxCombPC[currTF])]))
          }else{
            aggregPC.unscaled[ currTF , rownames(currWeightedPC) ] <-
              sqrt(sum(currWeightedPC[ , 1:(maxCombPC[currTF])]))
          }
          rm(currWeightedPC, currVarExp)
        }
        # Calculate the mean expression level of gene set
        mean_express<- list()
        for (currTF in names(tfTgtPca)){
          genes <- geneIds(predBind[[currTF]])
          mean_express[[currTF]] <- apply(mRNAExpr[which(rownames(mRNAExpr) %in% genes),],2,mean)
        }
        mean_express.1 <-mat.or.vec(nr = length(rownames(aggregPC.unscaled)), nc = length(colnames(aggregPC.unscaled)))
        rownames(mean_express.1) <- rownames(aggregPC.unscaled)
        colnames(mean_express.1) <- colnames(aggregPC.unscaled)
        for (currTF in rownames(aggregPC.unscaled)){
          mean_express.1[currTF,] <- mean_express[[currTF]]
        }
        ## ----scPS with scaled PC score above 0 and average expression of the gene set--------------------------------------------------
        new.spScore <- mat.or.vec(nr = length(rownames(aggregPC.unscaled)), nc = length(colnames(aggregPC.unscaled)))
        rownames(new.spScore) <- rownames(aggregPC.unscaled)
        colnames(new.spScore) <- colnames(aggregPC.unscaled)
        for (j in 1:length(colnames(mRNAExpr))){
          for (i in 1:length(names(tfTgtPca))){
            new.spScore[i,j] <- aggregPC.unscaled[i,j] * mean_express.1[i,j] 
          }
        }
        write.csv(new.spScore,file = paste("Score/rep",repl,'scPS.csv',sep = " "))
        }, 
      error = function(err) {
        print(paste("Error occurred in:", conditionMessage(err)))
        # cumulVarExp and CombPC
        varExplThresh <- 0.5 # target for cumulative max var explained
        cumulVarExp <- function(sdevPC){cumsum(sdevPC^2)/sum(sdevPC^2)}
        varExp <- function(sdevPC){(sdevPC^2) / sum(sdevPC^2)}
        maxCombPC <- unlist(lapply(tfTgtPca, FUN = function(x){
          idx = min(which(cumulVarExp(x[['pca']]@stdev) > varExplThresh));
          if(identical(idx, integer(0))){return(1)}else{return(max(idx))}
        }))
        table(maxCombPC)
        # Scale PC score above 0
        aggregPC.unscaled <- mat.or.vec(nr = length(names(tfTgtPca)),nc = ncol(mRNAExpr))
        rownames(aggregPC.unscaled) <- names(tfTgtPca)
        colnames(aggregPC.unscaled) <- colnames(mRNAExpr)
        PC.X <- list()
        for(currTF in rownames(aggregPC.unscaled)){
          PC.X[[currTF]] <- tfTgtPca[[currTF]]@reductions$pca@cell.embeddings - min(tfTgtPca[[currTF]]@reductions$pca@cell.embeddings)
        }
        # for each gene set after filtering  
        for(currTF in rownames(aggregPC.unscaled)){
          # get the var% explained per PC
          currVarExp = varExp(tfTgtPca[[currTF]]@reductions$pca@stdev)
          # weight PC scores
          currWeightedPC = sapply(1:maxCombPC[currTF], FUN =
                                    function(i){(PC.X[[currTF]][,i] * currVarExp[i])})
          if(ncol(currWeightedPC) > 1){
            aggregPC.unscaled[ currTF , rownames(currWeightedPC) ] <-
              sqrt(rowSums(currWeightedPC[ , 1:(maxCombPC[currTF])]))
          }else{
            aggregPC.unscaled[ currTF , rownames(currWeightedPC) ] <-
              sqrt(sum(currWeightedPC[ , 1:(maxCombPC[currTF])]))
          }
          rm(currWeightedPC, currVarExp)
        }
        # Calculate the mean expression level of gene set
        mean_express<- list()
        for (currTF in names(tfTgtPca)){
          genes <- geneIds(predBind[[currTF]])
          mean_express[[currTF]] <- apply(mRNAExpr[which(rownames(mRNAExpr) %in% genes),],2,mean)
        }
        mean_express.1 <-mat.or.vec(nr = length(rownames(aggregPC.unscaled)), nc = length(colnames(aggregPC.unscaled)))
        rownames(mean_express.1) <- rownames(aggregPC.unscaled)
        colnames(mean_express.1) <- colnames(aggregPC.unscaled)
        for (currTF in rownames(aggregPC.unscaled)){
          mean_express.1[currTF,] <- mean_express[[currTF]]
        }
        ## ----scPS with scaled PC score above 0 and average expression of the gene set--------------------------------------------------
        new.spScore <- mat.or.vec(nr = length(rownames(aggregPC.unscaled)), nc = length(colnames(aggregPC.unscaled)))
        rownames(new.spScore) <- rownames(aggregPC.unscaled)
        colnames(new.spScore) <- colnames(aggregPC.unscaled)
        for (j in 1:length(colnames(mRNAExpr))){
          for (i in 1:length(names(tfTgtPca))){
            new.spScore[i,j] <- aggregPC.unscaled[i,j] * mean_express.1[i,j] 
          }
        }
        write.csv(new.spScore,file = paste("Score/rep",repl,'scPS.csv',sep = " "))
      },
      finally = {
        message('All done, quitting.')
      }
    )
}
