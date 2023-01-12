library(getopt)
#library(Matrix)

arg <- matrix(c("input", "i","1","character","input file1",
                "outdir","o","1","character","outdir",
                "sample","s","1","character","sample,default=Maize",
                "tissue","t","1","character","tissue,default=Embro",
                "kfilt","k","1","integer","k.filter for merge",
                "dims","d","1","integer","dims option for FindNeighbors,default=15",
                "resolution","r","1","numeric","resolution option for FindClusters,[0.4-1.2],default=1",
                "help","h","0","logical", "Usage: Rscript runiDrop.R -i <input> -o <outdir> [-s SAMPLE -t TISSUE]",
                "minCG","m","1","integer","minimum gene number for each cell, i.e. nFeature_RNA, default=200",
                "rds","f","0","logical", "Save the RDS file"
               #"list","l","1","character", "interested gene"
                ),byrow=T,ncol=5)
opt = getopt(arg)
if(!is.null(opt$help) || is.null(opt$input)){
        cat(paste(getopt(arg, usage = T), "\n"))
 q()
}
if (is.null(opt$sample)){
        opt$sample <- "Maize"
}
if (is.null(opt$tissue)){
        opt$tissue <- "Embro"
}
if (is.null(opt$outdir)){
        opt$outdir <- "output"
}
if (is.null(opt$dims)){
        opt$dims <- 15
}
if (is.null(opt$resolution)){
        opt$resolution <- 1
}
if (is.null(opt$minCG)){
 opt$minCG <- 200
}


library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
mypalette<-c(brewer.pal(11,"BrBG"),brewer.pal(11,"PiYG"))

file=scan(opt$input,"")
dir.create(opt$outdir);
setwd(opt$outdir);
BRAIN<-list()
n=0
w=c()
h=c()
id=c()
reference<-readRDS("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhongliyuan/zhongliyuan/soybean_seed/soybean_root/testing_spatial_seurat_integrate/k50/data.rds")
reference<-FindVariableFeatures(reference)

for(i in file)
{
n=n+1
brain.data <- read.csv(i, header = T, sep=",")
rownames(brain.data) <- brain.data[,1]
brain.data <- brain.data[,2:dim(brain.data)[2]]
meta <- data.frame(Barcodes = colnames(brain.data), Sample = opt$sample, Tissue = opt$tissue, Annotation = NA, Celltype = NA)
BRAIN[[n]] <- CreateSeuratObject(counts = brain.data,
                              project = opt$tissue,
                              min.cells = 0,
                              min.features = 0,
                              assay = 'Spatial')

#adding location information
locations<-data.frame(substring(colnames(brain.data),2))
locations[,1:2]<-str_split_fixed(locations$substring.colnames.brain.data...2.,"_",2)
names(locations)<-c("x","y")
rownames(x = locations) <-colnames(brain.data)
locations$x<-as.numeric(locations$x)
locations$y<-as.numeric(locations$y)
w[n]<-as.numeric((max(locations$x) - min(locations$x) + 1)/10)
h[n]<-as.numeric((max(locations$y) - min(locations$y) + 1)/10)
BRAIN[[n]][['image']] <- new(Class = 'SlideSeq',assay = "Spatial",
    coordinates = locations)
    name<-strsplit(i,"/",fixed=T)[[1]][2]
    BRAIN[[n]]@meta.data$group<-strsplit(name,".",fixed=T)[[1]][1]
    id[n]<-strsplit(name,".",fixed=T)[[1]][1]
    BRAIN[[n]] <- SCTransform(BRAIN[[n]], assay = "Spatial", verbose = FALSE) 
    BRAIN[[n]] <- FindVariableFeatures(BRAIN[[n]], selection.method = "vst", nfeatures = 3000) 
    BRAIN[[n]] <- RunPCA(object = BRAIN[[n]], verbose = FALSE)
    BRAIN[[n]] <- RunUMAP(BRAIN[[n]], dims = 1:opt$dims)
 
    anchors <- FindTransferAnchors(
               reference = reference,
               query = BRAIN[[n]],
               dims=1:opt$dims)
    predictions <- TransferData(
                   anchorset = anchors,
                   refdata = Idents(reference),
                   dims=1:opt$dims
                   )
    BRAIN[[n]] <- AddMetaData(BRAIN[[n]], metadata = predictions)
    pdf(paste0(BRAIN[[n]]@meta.data$group,"_UMAP.pdf"), width = 10)
    plot1 <- DimPlot(object = BRAIN[[n]], reduction = "umap", group.by="predicted.id",label = T,pt.size=1)
    print(plot1)
    dev.off()

    pdf(paste0(BRAIN[[n]]@meta.data$group,"_Dimspatial.pdf"))
    plot2 <- SpatialDimPlot(BRAIN[[n]], pt.size=1,group.by="predicted.id",label = T, label.size = 1)+scale_fill_manual(breaks=levels(Idents(reference)),values=mypalette)
    print(plot2)
    dev.off()
    saveRDS(BRAIN[[n]],file=paste0(BRAIN[[n]]@meta.data$group,"_Data.RDS")[1])
}














