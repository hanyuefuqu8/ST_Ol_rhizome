library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(readr)
library(cowplot)
BRAIN<-readRDS("../monocle3_2/sub_data_merge.RDS")
rdsf <- read_tsv("../info3")

temp<-left_join(BRAIN@meta.data,rdsf,by="group")
BRAIN@meta.data$phase<-temp$phase
BRAIN@meta.data$cell.type<-paste0("phase",BRAIN@meta.data$phase,"_",Idents(BRAIN))

#change the first name of images
b<-paste0("image_",names(table(BRAIN$group)))
c<- names(BRAIN@images)
names(BRAIN@images)[1]<-b[!b%in%c]

marker1<-FindMarkers(BRAIN,
                     ident.1="phase1_9",
					 ident.2="phaseck_2",
					 group.by="cell.type")
write.table(file="p1_c9_ck2.txt",marker1,sep="\t")

Idents(BRAIN)<-BRAIN$cell.type
gene1 <- rownames(marker1[marker1$avg_log2FC>0 & marker1$p_val<0.01,])

marker2<-read.table("../Diff_genes_2/p1_phase1_3_phaseck_0.txt")
gene2 <- rownames(marker2[marker2$avg_log2FC>0 & marker2$p_val<0.01,])

#bottom <- marker1  %>% top_n(n = -5, wt = avg_log2FC)
gene<-c(gene1,gene2)
pdf("p1_c9_c3_ck2_ck0_dotplot.pdf",height=5)
DotPlot(BRAIN, features =unique(gene),assay="SCT",group.by="cell.type",idents=c("phase1_9","phase1_3","phaseck_2","phaseck_0"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

for (i in gene)
{
library(gridExtra)
p<-list()
for(j in names(BRAIN@images))
{
p[[j]]<- SpatialFeaturePlot(BRAIN,images=j,features=i,stroke=0,pt.size=12)+ggtitle(j)
}
sum<-grid.arrange(grobs=p,ncol=4)
ggsave2(paste0(i,"sub_dimspatial.pdf"),sum,width=20,height=20)

}

