library(Seurat)
library(SeuratDisk)
library(ggplot2)

file<-scan("file_and_ck","")
#取出芽处的所有细胞
cell<-c()
for (i in file)
{
Convert(i,assay='SCT',dest = "h5seurat", overwrite = TRUE)
name<-paste0(strsplit(i,'.',fixed=T)[[1]][1],".h5seurat")
data <- LoadH5Seurat(name)
if(i=="BIN100165TR_A2down02ck.h5ad")
{
cells<-paste0(strsplit(i,'ck.',fixed=T)[[1]][1],"_X",data@meta.data$x,"_",data@meta.data$y)
cellck<-cells
} else
{
cells<-paste0(strsplit(i,'.',fixed=T)[[1]][1],"_X",data@meta.data$x,"_",data@meta.data$y)
}
cell<-c(cell,cells)
}

#取子集
sdata<-readRDS("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/baiyixuan/Ol_data/baiyixuan/new_stage/harmony/data.RDS")
BRAIN<-subset(sdata,cells=cell)
BRAIN@meta.data[rownames(BRAIN@meta.data)%in%cellck,'group']="BIN100165TR_A2down02ck"
BRAIN@images<-BRAIN@images[-4]
saveRDS(BRAIN,"sub_data_and_ck.RDS")
