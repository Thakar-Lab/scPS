setwd("/Users/rwang4/Library/CloudStorage/Box-Box/Ruoqiao Thakar Lab/Single-cell Pathway Score (scPS)/Data & results/HIV and B cell/")
if(!require(Seurat))install.packages("Seurat")
if(!require(ggplot2))install.packages("ggplot2")
if(!require(dplyr))install.packages("dplyr")
library(ggsci)
library(ggrepel)
rm(list=ls()) 
gc()

My_Theme <- theme(axis.ticks = element_blank(),
                  legend.title=element_blank(),
                  legend.position = "bottom", 
                  legend.direction = "horizontal", 
                  legend.justification='center',
                  plot.caption = element_text(hjust = 0.5, vjust = -1, size=14))
# Started here 
###########################################################################################################################################
# TSNE plot
# load HIV_data.RData
load("~/Library/CloudStorage/Box-Box/Ruoqiao Thakar Lab/Single-cell Pathway Score (scPS)/Data & results/HIV and B cell/HIV_data.RData")
Idents(Participant.seurat) <- "ClusterName"
DimPlot(Participant.seurat, reduction = "tsne",label = T, pt.size = 0.8) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 10), legend.text = element_text(size = 10)) + 
  guides(colour = guide_legend(override.aes = list(size = 5))) +NoLegend()

umap = Participant.seurat@reductions$tsne@cell.embeddings %>%
  as.data.frame() %>% 
  cbind(cell_type = Participant.seurat@meta.data$ClusterName)

cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2))

# set the colour
colour <- pal_d3("category20",alpha = 0.6)(16)
cluster <- unique(umap$cell_type)
cluster <- cluster[order(cluster)]
names(colour) <- cluster

p <- ggplot(umap,aes(x= tSNE_1 , y = tSNE_2 ,color = cell_type)) +
  geom_point(size = 1 , alpha =0.8 )  +
  scale_color_manual(values = colour)+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_rect(colour = "black",fill = NA), #边框
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=10), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小hiv
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中点的大小
  geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines"))+
  theme(legend.position = "none")
p
# pdf("tsne_cluster.pdf",height = 10, width = 10)
# p
# dev.off()

library(grid)
library(gridExtra)
library(ggpubr)

new.spScore <- read.csv("scPS for HIV-data.csv",header = T,row.names = 1)
AUCell <- read.csv("AUCell for HIV data.csv",header = T,row.names = 1)
JASMINE <- read.csv("JASMINE for HIV data.csv",header = T,row.names = 1)
UCell <- read.csv("UCell for HIV data.csv",header = T,row.names = 1)
ssGSEA <- read.csv("ssGSEA for HIV data.csv",header = T,row.names = 1)
SCSE<- read.csv("scSE for HIV data.csv",header = T,row.names = 1)
AddModuleScore<- read.csv("AddModuleScore for HIV data.csv",header = T,row.names = 1)

my_data = SCSE
#######################################################################################
# Apply the scale function to each row from row 1 to row 50 
scale_row_to_01 <- function(row) {
  (row - min(row, na.rm = TRUE)) / (max(row, na.rm = TRUE) - min(row, na.rm = TRUE))
}
# Apply the scale function to each row from row 1 to row 50
scaled_data <- my_data
for (i in 1:dim(my_data)[1]) {
  scaled_data[i,] <- scale_row_to_01(my_data[i,])
}
dat <- cbind(umap,t(scaled_data))
#
######### select all the cell type #########################
# from above
cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2))

######### select the cell type #########################
# dat$cell_type <- gsub("B cells naive - 1|B cells naive - 2|B cells memory", "B cells", dat$cell_type)
# dat$cell_type <- gsub("NK cells resting|NK cells activated", "NK cells", dat$cell_type)
# exclude_cell_type <- c("B cells","NK cells","Monocytes")
# dat$cell_type <- gsub(paste0("^(?!",paste(exclude_cell_type,collapse = "|"),"$).*$"),"T cells", dat$cell_type,perl=T)
#
# cell_type_med$cell_type  <- gsub("B cells naive - 1|B cells naive - 2|B cells memory", "B cells", cell_type_med$cell_type )
# cell_type_med$cell_type  <- gsub("NK cells resting|NK cells activated", "NK cells", cell_type_med$cell_type )
# cell_type_med$cell_type[7:16]  <- "T cells"
# cell_type_med<- cell_type_med[c(3,4,6,8),]

dat$cell_type <- gsub("B cells naive - 1|B cells naive - 2", "Naive B", dat$cell_type)
dat$cell_type <- gsub("B cells memory", "Memory B", dat$cell_type)
dat$cell_type <- gsub("T cells CD8 - 1|T cells CD8 - 2|T cells CD8 - 3|T cells CD8 - 4", "CD8 T", dat$cell_type)
dat$cell_type <- gsub("T cells CD4 - 1", "CD4 T", dat$cell_type)
dat$cell_type <- gsub("T cells naive", "Naive T", dat$cell_type)
dat$cell_type <- gsub("T cells CD8/NK cells resting", "CD8 T/resting NK", dat$cell_type)
dat$cell_type <- gsub("T cells CD8/CD4/CD4 naive", "CD8/CD4/naive T", dat$cell_type)
dat$cell_type <- gsub("T cells CD8 naive/T cells CD8/ RBCs", "CD8/naive T/RBCs", dat$cell_type)
dat$cell_type <- gsub("T cells CD4 memory resting/ T cells CD8", "CD8/CD4 memory T", dat$cell_type)
dat$cell_type <- gsub("NK cells activated", "Activated NK", dat$cell_type)
dat$cell_type <- gsub("NK cells resting", "Resting NK", dat$cell_type)
unique(dat$cell_type)

cell_type_med$cell_type <- gsub("B cells naive - 1|B cells naive - 2", "Naive B", cell_type_med$cell_type)
cell_type_med$cell_type <- gsub("B cells memory", "Memory B", cell_type_med$cell_type)
cell_type_med$cell_type <- gsub("T cells CD8 - 1|T cells CD8 - 2|T cells CD8 - 3|T cells CD8 - 4", "CD8 T", cell_type_med$cell_type)
cell_type_med$cell_type <- gsub("T cells CD4 - 1", "CD4 T", cell_type_med$cell_type)
cell_type_med$cell_type <- gsub("T cells naive", "Naive T", cell_type_med$cell_type)
cell_type_med$cell_type <- gsub("T cells CD8/NK cells resting", "CD8 T/resting NK", cell_type_med$cell_type)
cell_type_med$cell_type <- gsub("T cells CD8/CD4/CD4 naive", "CD8/CD4/naive T", cell_type_med$cell_type)
cell_type_med$cell_type <- gsub("T cells CD8 naive/T cells CD8/ RBCs", "CD8/naive T/RBCs", cell_type_med$cell_type)
cell_type_med$cell_type <- gsub("T cells CD4 memory resting/ T cells CD8", "CD8/CD4 memory T", cell_type_med$cell_type)
cell_type_med$cell_type <- gsub("NK cells activated", "Activated NK", cell_type_med$cell_type)
cell_type_med$cell_type <- gsub("NK cells resting", "Resting NK", cell_type_med$cell_type)

#### reorder the cell type ####
unique(dat$cell_type)
dat$cell_type <- factor(dat$cell_type, levels = c("CD8 T",
                                                  "CD4 T",
                                                  "Naive T",
                                                  "CD8/CD4/naive T",
                                                  "CD8/CD4 memory T",
                                                  "CD8 T/resting NK",
                                                  "CD8/naive T/RBCs",
                                                  "Naive B",
                                                  "Memory B",
                                                  "Monocytes",
                                                  "Activated NK",
                                                  "Resting NK"))

cell_type_med$cell_type <- factor(cell_type_med$cell_type, levels = c("CD8 T",
                                                                      "CD4 T",
                                                                      "Naive T",
                                                                      "CD8/CD4/naive T",
                                                                      "CD8/CD4 memory T",
                                                                      "CD8 T/resting NK",
                                                                      "CD8/naive T/RBCs",
                                                                      "Naive B",
                                                                      "Memory B",
                                                                      "Monocytes",
                                                                      "Activated NK",
                                                                      "Resting NK"))

#########  select the GS ##############################
GS_score <-c()
GS <- colnames(dat)[4:16]
GS <- GS[c(3,8,10,11)]
#######################################################
# boxplot
# for all cell type
# set the colour
colour <- pal_d3("category20",alpha = 1)(16)
cluster <- unique(dat$cell_type)
cluster <- cluster[order(cluster)]
names(colour) <- cluster

currGS=GS[4]
currGS
unique(dat$cell_type)
c("Naive B","CD8 T","Monocytes","Resting NK")
for (currGS in GS) {
  # save 6x8 inches, landscape
  p3 <- ggplot(dat, aes(x=cell_type, y= get(currGS), color=cell_type)) +
    geom_boxplot()+
    stat_compare_means(label="p.signif",ref.group = 'CD8 T') +     
    scale_color_manual(values=colour)+
    theme_bw(base_size = 14)+ 
    labs(x ="", y = paste("SCSE of T gene set")) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "right",
          axis.ticks = element_blank(),
          legend.title=element_blank(),
          plot.caption = element_text(hjust = 0.5, vjust = -1, size=10))+
    theme(legend.position = "none") 
  
  p3
  
  ggsave(p3, file=paste0("/Users/rwang4/Desktop/",deparse(substitute(SCSE)),"_",currGS,".pdf"), width =8, height = 6)
 }

# # for selected cell type
# for (currGS in GS) {
#   p3 <- ggplot(dat, aes(x=cell_type, y= get(currGS), color=cell_type)) +
#   geom_boxplot()+
#   #stat_compare_means(ref.group = GS_cell_type[which(GS ==currGS)]) +                     
#   scale_color_manual(values=colour)+
#   theme_bw(base_size = 14)+ 
#   labs(x ="", y = paste("scPS of",GS_score[which(GS ==currGS)])) + 
#     theme(axis.text.x = element_blank(),
#         legend.position="right",
#         axis.ticks = element_blank(),
#         legend.title=element_blank(),
#         plot.caption = element_text(hjust = 0.5, vjust = -1, size=14))
#   
#   
#   p3
#   
#   ggsave(p, file=paste0("./Figures/Boxplot_",deparse(substitute(ssGSEA)),"_",GS_cell_type[which(GS ==currGS)],".pdf"), width = 6, height = 4)
# }

#### Feature plot ####
currGS=GS[4]
currGS
for (currGS in GS) {
  p4 = ggplot(dat, aes(x = tSNE_1, y = tSNE_2, color = get(currGS))) +
    geom_point(size = 1, alpha = 1) +
    scale_color_gradient(low = "white", high = "darkblue") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = 'white'),
      plot.background = element_rect(fill = "white")
    ) +
    theme(
      legend.title = element_blank(),
      legend.key = element_rect(fill = 'white'),
      legend.text = element_text(size = 16),
      legend.key.size = unit(1, 'cm')
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))+  
    geom_label(
      data = cell_type_med,
      aes(x = tSNE_1, y = tSNE_2, label =cell_type),
      size=5,
      color = "black")+
    theme(legend.position = "none") 
  

  p4
  
  #ggtitle(currGS)
  ggsave(p4, file=paste0("/Users/rwang4/Desktop/",deparse(substitute(SCSE)),"_feature_",currGS,".pdf"), width = 10, height = 8)
}


p <- ggplot(dat,aes(x= tSNE_1 , y = tSNE_2 ,color = cell_type)) +
  geom_point(size = 1 , alpha =0.5 )  +
  scale_color_manual(values = colour)+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_rect(colour = "black",fill = NA), #边框
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))+
  theme(legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=10), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小hiv
  guides(color = guide_legend(override.aes = list(size=5)))+ #设置legend中点的大小
  geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines"))+
  theme(legend.position = "none")
p















































































































































