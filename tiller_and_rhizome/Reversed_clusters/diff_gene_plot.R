library(getopt)
#library(Matrix)

arg <- matrix(c("input", "i","1","character","input file1",
		"geneid", "g","1","character","table include geneid and genename",
		"info", "f","1","character","information file",
		"method","m","1","character","method to normalize data"
                ),byrow=T,ncol=5)
opt = getopt(arg)

library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(cowplot)
if(FALSE)
{
opt<-list()
opt$input<-"../data.RDS"
opt$geneid<-"geneid"
opt$info<-"../diff_gene_SCT/info"
opt$method<-"SCT"
}
BRAIN<-readRDS(opt$input)
genes<- read.table(opt$geneid,header=T)

if(opt$method=="SCT")
{
#BRAIN.list <- SplitObject(BRAIN, split.by = "group")
#for( i in 1:length(BRAIN.list))
#{
#BRAIN.list[[i]] <- SCTransform(BRAIN.list[[i]], vst.flavor = "v2",assay='Spatial',verbose = FALSE)
#}
#BRAIN<-merge(BRAIN.list[[1]], y=BRAIN.list[-1])
#BRAIN<-PrepSCTFindMarkers(BRAIN)
DefaultAssay(BRAIN) <- "SCT"} else {
DefaultAssay(BRAIN) <- "Spatial"
BRAIN <- NormalizeData(BRAIN, normalization.method = "LogNormalize", scale.factor = 10000)
}

names(BRAIN@images)

dir.create("images")
setwd("images")


for(j in genes[,1])
{
library(gridExtra)
p<-list()
exp<-data.frame(FetchData(object = BRAIN, vars = j))
n=1
for(i in names(table(BRAIN$group)))
{
name<-gsub("-",".",i)
#if(i!=names(table(BRAIN$group))[1])
#{
#imagename<-paste0("image_",i,".",n)
#n=n+1} else {
#imagename<-paste0("image_",i)
#}

imagename<-paste0("image_",name)

temp<-merge(BRAIN@images[[imagename]]@coordinates[,1:2],exp,by = 'row.names')
colnames(temp)[4]<-"expression"
#h <- as.numeric(max(temp$y) - min(temp$y) + 1)
#w <- as.numeric(max(temp$x) - min(temp$x) + 1)
p[[i]]<-ggplot(temp) + geom_point(mapping=aes(x=x, y=y, colour = expression),size = 6)+
  scale_color_gradient(limit=c(0,max(exp)),low = "blue",high = "#eaed18")+
  ggtitle(i)

}
sum<-grid.arrange(grobs=p,ncol=4)
ggsave2(paste0(j,".pdf"),sum,width=25,height=20)

}

setwd("../")
rdsf <- read.table(opt$info,sep="\t",header=T)

temp<-left_join(BRAIN@meta.data,rdsf,by="group")
BRAIN@meta.data$chip_name<-temp$chip_name
BRAIN@meta.data$line<-temp$line
BRAIN@meta.data$species<-temp$species
BRAIN@meta.data$cell.type<-Idents(BRAIN)
BRAIN$celltype.treat <- paste(BRAIN$species, Idents(BRAIN), sep = "_")

temp<-genes[!duplicated(genes[,1]),]
temp[,2]<-make.unique(temp[,2],sep="_")

if(opt$method=="SCT")
{a<-temp[,1]%in%rownames(BRAIN[["SCT"]]@data)} else 
{a<-temp[,1]%in%rownames(BRAIN[["Spatial"]]@data)}

temp<-temp[a,]
#dotplot
pdf("diff_gene_dot_plot.pdf",height=10,width=12)
p1<-DotPlot(BRAIN, features = unique(temp[,1]),cols = c("blue", "red"), dot.scale = 8, split.by = "species")+scale_x_discrete(labels=temp[,2])+theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1
dev.off()



#Vlnplot
pdf("diff_gene_Vln_plot.pdf",height=5,width=10)
p1<-VlnPlot(BRAIN, features = unique(temp[,1]),group.by="cell.type",split.by = "species",pt.size = 0, combine = FALSE)
p1
dev.off()


#heatmap
Idents(BRAIN) <- "celltype.treat"
BRAIN <- ScaleData(BRAIN)
pdf("diff_gene_heatmap_plot.pdf", width = 6, height = 10)
BRAINs<-subset(x=BRAIN,idents=c("Aerial_stem_4","Aerial_stem_5","Underground_stem_4","Underground_stem_5"))
DoHeatmap(BRAINs,features = unique(temp[,1]))+scale_fill_gradientn(colors = c( "red","black", "green"))
dev.off()


p<-AverageExpression(BRAIN,features=temp[,1],slot='counts',group.by = "ident")

if (opt$method=="SCT")
{
temp3<-as.data.frame(p$SCT)
}else 
{
temp3<-as.data.frame(p$Spatial)
}

temp3$geneid<-rownames(temp3)
temp4<-merge(temp3,temp,by="geneid")
rownames(temp4)<-temp4$genename
num=ncol(temp4)
temp4<-temp4[,-c(1,num)]
temp4<-temp4[,order(colnames(temp4))]


p1<-pheatmap::pheatmap(temp4,
                       scale="row",
                       cluster_cols= FALSE,
                       cluster_rows=FALSE,
                       #gaps_col = seq(from=2,by=2,length=8),
                       )
ggsave2("diffgene_heatmap.pdf",p1,height=14,width=8)

