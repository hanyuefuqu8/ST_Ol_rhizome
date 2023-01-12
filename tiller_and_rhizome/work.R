library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
#library(gdata)
pancreas.integrated<-readRDS("data.RDS")

if(TRUE)
{
pancreas.integrated<-RenameIdents(pancreas.integrated, 
`0` = "Vascular bundle 1",
 `1` = "Parenchyma 1",
 `2` = "Leaf promodium 1",
 `3` = "Leaf promodium 2",
 `4` = "Vascular bundle 2",
 `5` = "Parenchyma 2",
 `6` = "Parenchyma 3",
 `7` = "Meristem 1",
 `8` = "Apical initial cells",
 `9` = "Unknown",
`10` = "Sclerenchyma",
`11` = "Epidermis",
`12` = "Meristem 2",
`13` = "Vascular bundle 3",
`14` = "Aerenchyma",
`15` = "Vascular bundle 4",
`16` = "Apical meristem",
`17` = "Unknown"
)
#markerlist<-read.table("markergene.lst",sep="\t",header=T)

info<-read.table("markerlist",header=T,sep="\t")
temp<-info[!duplicated(info$geneid),]
a<-temp[,'geneid']%in%rownames(pancreas.integrated[["SCT"]]@data)
temp<-temp[a,]
markerlist<-temp

}

#reset colors
library(gridExtra)
library(RColorBrewer)
#mypalette<-c("#5cb3cc","#1a6840","#ef6f48","#0f95b0","#2775b6","#1a94bc","#bacf65","#525288","#d276a3","#c08eaf","#fffef9")
mypalette<-c("#0f95b0","#9fa39a","#5cb3cc","#bacf65","#d276a3","#e2e7bf","#add5a2","#5dbe8a","#229453","#253d24","#1a6840","#f17666","#525288","#c08eaf","#d2b116","#f59f89","#f7de98")
#change the first name of images
#b<-paste0("image_",names(table(pancreas.integrated$group)))
#c<- names(pancreas.integrated@images)
#names(pancreas.integrated@images)[1]<-b[!b%in%c]
#展示抠出来的图
library(gridExtra)
p<-list()
for(i in names(pancreas.integrated@images))
{
p[[i]]<- SpatialDimPlot(pancreas.integrated,images=i,stroke=0,pt.size=12,label=F)+scale_fill_manual(breaks=levels(Idents(pancreas.integrated)),values=mypalette)+ggtitle(i)
}
sum<-grid.arrange(grobs=p,ncol=4)
ggsave2("dimspatial.pdf",sum,width=20,height=20)

pdf("09_umap_3.pdf", width = 10,height=8)
plot1 <- DimPlot(object = pancreas.integrated, reduction = "umap", label = T,label.size=6,pt.size=1,group.by = "ident",cols=mypalette)
plot1
dev.off()

file=scan("file","")
Idents(pancreas.integrated)<- ordered(Idents(pancreas.integrated), levels =c( "Aerenchyma",
 "Apical initial cells",
 "Apical meristem",
 "Meristem 1",
 "Meristem 2",
 "Parenchyma 2",
 "Parenchyma 1",
 "Parenchyma 3",
 "Leaf promodium 1",
 "Leaf promodium 2",
 "Sclerenchyma",
 "Vascular bundle 3",
 "Vascular bundle 4",
 "Vascular bundle 1",
 "Vascular bundle 2",
 "Epidermis",
 "Unknown"
)
 )

p2<-DotPlot(pancreas.integrated, features = markerlist[,'geneid'],assay="SCT")+
            scale_x_discrete(labels=markerlist[,'genename'])+
			theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y=element_text(face="italic"))+
			scale_colour_gradient(low = "#fa87bf", high = "#1c07a6")+
                        coord_flip()
ggsave2("selected_marker_bubble.pdf",p2,dpi=300,width=6,height=6)

#
#BRAIN.list <- SplitObject(pancreas.integrated, split.by = "group")
#for(i in 2:length(file))
#{
 #  BRAIN.list[[i]]@images<- BRAIN.list[[i]]@images[i]
  # pdf(paste0(BRAIN.list[[i]]@meta.data$group,"_dimspatial.pdf"))
   #plot1 <- SpatialDimPlot(BRAIN.list[[i]],stroke=0,pt.size=8,label=T,label.size=3)+scale_fill_manual(breaks=levels(Idents(pancreas.integrated)),values=mypalette)
   #print(plot1)
   #dev.off()

   #pdf(paste0(BRAIN.list[[i]]@meta.data$group,"_dimsplit.pdf"), width = 10,height=10)
   #plot1 <-SpatialDimPlot(BRAIN.list[[i]], cells.highlight = CellsByIdentities(object = BRAIN.list[[i]], idents = levels(Idents(BRAIN.list[[i]]))), facet.highlight = TRUE, ncol = 5, stroke = 0, pt.size=1)
   #print(plot1)
   #dev.off() 
#}



saveRDS(pancreas.integrated,"data2.RDS")




