library(Seurat)
BRAIN<-readRDS("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/yangfeng/baiyixuan/new_stage/harmony/data.RDS")
BRAINi<-SplitObject(BRAIN,split.by="group")

for(i in 1:length(BRAINi))
{

temp<-cbind(BRAINi[[i]]@images[[i]]@coordinates[,1:2],as.numeric(BRAINi[[i]]@meta.data[,'seurat_clusters'])*100)

temp$geneID<-rownames(temp)
colnames(temp)[3]<-"MIDCounts"
 write.table(file=paste0(BRAINi[[i]]$group[1],".gem"),temp[,c('geneID','x','y','MIDCounts')],row.names=F,quote=F)
 }
