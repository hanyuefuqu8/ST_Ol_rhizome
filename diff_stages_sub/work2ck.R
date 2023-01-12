library(Seurat)
library(SeuratDisk)
library(ggplot2)

file<-scan("fileck","")
#取出芽处的所有细胞
cell<-c()
for (i in file)
{
Convert(i, dest = "h5seurat", overwrite = TRUE)
name<-paste0(strsplit(i,'.',fixed=T)[[1]][1],".h5seurat")
data <- LoadH5Seurat(name)
cells<-paste0(strsplit(i,'ck.',fixed=T)[[1]][1],"_X",data@meta.data$x,"_",data@meta.data$y)
cell<-c(cell,cells)
}

#取子集
sdata<-readRDS("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/yangfeng/baiyixuan/new_stage/harmony/data.RDS")
BRAIN<-subset(sdata,cells=cell)
BRAIN@images<-BRAIN@images[-4]
saveRDS(BRAIN,"sub_datack.RDS")
