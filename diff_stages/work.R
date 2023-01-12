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
pancreas.integrated<-RenameIdents(pancreas.integrated, `0` = "Parenchyma 1", `1` = "Sclerenchyma", `2` = "Parenchyma 2",`3` = "Leaf promeristem", `4` = "Phloem", `5` = "Xylem parenchyma", `6` = "Epidermis 1", `7` = "Aerenchyma", `8` = "Unknown", `9` = "Meristem",`10` = "Epidermis 2")
Idents(pancreas.integrated)<- ordered(Idents(pancreas.integrated), levels =c(
 "Parenchyma 1",
 "Sclerenchyma",
 "Parenchyma 2",
 "Meristem",
 "Leaf promeristem",
 "Phloem",
 "Xylem parenchyma",
 "Epidermis 1",
 "Epidermis 2",
 "Aerenchyma",
 "Unknown"
) )
#markerlist<-read.table("markergene.lst",sep="\t",header=T)

info<-read.table("markerlist",header=T,sep="\t")
temp<-info[!duplicated(info$geneid),]
a<-temp[,'geneid']%in%rownames(pancreas.integrated[["SCT"]]@data)
temp<-temp[a,]
markerlist<-temp

p2<-DotPlot(pancreas.integrated, features = markerlist[,'geneid'],assay="SCT")+
            scale_x_discrete(labels=markerlist[,'genename'])+
			theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.y=element_text(face="italic"))+
			scale_colour_gradient(low = "#fa87bf", high = "#1c07a6")+
                        coord_flip()
ggsave2("selected_marker_bubble.pdf",p2,dpi=300,width=6,height=6)
}

#reset colors
library(gridExtra)
library(RColorBrewer)
#mypalette<-c("#5cb3cc","#1a6840","#ef6f48","#0f95b0","#2775b6","#1a94bc","#bacf65","#525288","#d276a3","#c08eaf","#fffef9")
mypalette<-c("#5cb3cc","#2775b6","#0f95b0","#d276a3","#1a6840","#1a94bc","#d7e637","#ef6f48","#2775b6","#c08eaf","#f4f788")
#change the first name of images
#b<-paste0("image_",names(table(pancreas.integrated$group)))
#c<- names(pancreas.integrated@images)
#names(pancreas.integrated@images)[1]<-b[!b%in%c]
#展示抠出来的图
library(gridExtra)
p<-list()
for(i in names(pancreas.integrated@images))
{
p[[i]]<- SpatialDimPlot(pancreas.integrated,images=i,stroke=0,pt.size=6,label=F)+scale_fill_manual(breaks=levels(Idents(pancreas.integrated)),values=mypalette)+ggtitle(i)
}
sum<-grid.arrange(grobs=p,ncol=4)
ggsave2("dimspatial.pdf",sum,width=20,height=20)

pdf("09_umap_3.pdf", width = 12,height=8)
plot1 <- DimPlot(object = pancreas.integrated, reduction = "umap", label = T,label.size=7,pt.size=1,group.by = "ident",cols=mypalette)
plot1
dev.off()

file=scan("file","")

BRAIN.list <- SplitObject(pancreas.integrated, split.by = "group")
for(i in 2:length(file))
{
   BRAIN.list[[i]]@images<- BRAIN.list[[i]]@images[i]
   pdf(paste0(BRAIN.list[[i]]@meta.data$group,"_dimspatial.pdf"))
   plot1 <- SpatialDimPlot(BRAIN.list[[i]],stroke=0,pt.size=8,label=T,label.size=3)+scale_fill_manual(breaks=levels(Idents(pancreas.integrated)),values=mypalette)
   print(plot1)
   dev.off()

   pdf(paste0(BRAIN.list[[i]]@meta.data$group,"_dimsplit.pdf"), width = 10,height=10)
   plot1 <-SpatialDimPlot(BRAIN.list[[i]], cells.highlight = CellsByIdentities(object = BRAIN.list[[i]], idents = levels(Idents(BRAIN.list[[i]]))), facet.highlight = TRUE, ncol = 5, stroke = 0, pt.size=1)
   print(plot1)
   dev.off() 
}



#saveRDS(pancreas.integrated,"data2.RDS")




