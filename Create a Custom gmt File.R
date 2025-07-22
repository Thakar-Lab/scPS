# Create your own gmt file
# Create a Custom .gmt File from Harmonizome
# https://maayanlab.cloud/Harmonizome/
library(jsonlite)

get_gene_set <- function(term, dataset = "GeneRIF+Biological+Term+Annotations") {
  base_url <- "https://maayanlab.cloud/Harmonizome/api/1.0/gene_set"
  url <- paste0(base_url, "/", term, "/", dataset)
  response <- fromJSON(url)
  genes <- response[["associations"]]$gene$symbol
  return(genes)
}

macrophage_genes <- get_gene_set("macrophage")
monocyte_genes   <- get_gene_set("monocyte")

gene_sets <- list(
  macrophage = macrophage_genes,
  monocyte   = monocyte_genes
)

file = "Harmonizome_Monocyte_Macrophage.gmt"
write.gmt <- function(gs,file){
  sink(file)
  lapply(names(gs),function(i){
    cat(paste(c(i,'tmp',gs[[i]]),collapse = '\t'))
    cat('\n')
  })
  sink()
}
write.gmt(gene_sets,file)

# Check the gmt file
library(GSEABase)
library(genefilter)
predBind <- GSEABase::getGmt("Harmonizome_Monocyte_Macrophage.gmt")
print(length(predBind))