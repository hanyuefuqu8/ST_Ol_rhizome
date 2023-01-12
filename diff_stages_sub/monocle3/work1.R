library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)


BRAIN<-list()
BRAIN[[1]]<-readRDS("../sub_data.RDS")
BRAIN[[2]]<-readRDS("../sub_datack.RDS")
BRAIN[[2]]@meta.data$group<-"BIN100165TR_A2down02ck"

temp<-list()
temp[['image_BIN100165TR_A2down02ck']]<-BRAIN[[2]]@images[['image_BIN100165TR_A2down02']]
BRAIN[[2]]@images<-temp

pancreas.integrated<-merge(BRAIN[[1]], y=BRAIN[-1])
saveRDS(pancreas.integrated,"sub_data_merge.RDS")

BRAINs<-subset(x=pancreas.integrated,idents=c("0","4","2","3","5","9"))
saveRDS(BRAINs,"sub_sub_data_merge.RDS")

