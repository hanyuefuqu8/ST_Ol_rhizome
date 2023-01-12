

library(getopt)
#library(Matrix)

arg <- matrix(c("input", "i","1","character","input file1",
                "info", "f","1","character","information file",
		        "node1", "n","1","character","cluster contains node",
                        "node2", "m","1","character","cluster contains node",
		        "graph","p","1","character","choose knn or principal_graph"
                ),byrow=T,ncol=5)
opt = getopt(arg)



library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(monocle3)
library(readr)
library(cowplot)
library(viridis)
if(FALSE)
{
opt<-list()
opt$input<-"sub_sub_data_merge.RDS"
opt$info<-"../info3"
opt$node1<-"phaseck_Parenchyma"
opt$node2<-"phaseck_Xylem_parenchyma"
opt$graph<-"principal_graph"
}

BRAIN<-readRDS(opt$input)
BRAIN<-RenameIdents(BRAIN, `0` = "Parenchyma", `1` = "Specialized_parenchyma", `2` = "Peripheral_cylinder_of_vascular_bundles",`3` = "Leaf_promeristem", `4` = "Phloem", `5` = "Xylem_parenchyma", `6` = "Cuticle", `7` = "Aerenchyma", `8` = "Unknown", `9` = "Promeristem",`10` = "Epidermis")
Idents(BRAIN)<- ordered(Idents(BRAIN), levels =c("Parenchyma","Specialized_parenchyma","Peripheral_cylinder_of_vascular_bundles","Phloem", "Xylem_parenchyma","Aerenchyma","Leaf_promeristem","Promeristem","Cuticle","Epidermis", "Unknown") )
mypalette<-c("#5cb3cc","#0f95b0","#d276a3","#1a6840","#d7e637","#ef6f48")
rdsf <- read_tsv(opt$info)
temp<-left_join(BRAIN@meta.data,rdsf,by="group")
BRAIN@meta.data$phase<-temp$phase
BRAIN@meta.data$cell.type<-paste0("phase",BRAIN@meta.data$phase,"_",Idents(BRAIN))

#subset1
#BRAINs1<-subset(x=BRAIN,idents=c("Procambium","Xylem","Preprocambium","Stele initials"))
data <- as(as.matrix(BRAIN@assays$Spatial@counts), 'sparseMatrix')
pd <- data.frame(BRAIN@meta.data,Idents(BRAIN))
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(data,
                         cell_metadata  = pd,
                         gene_metadata  = fData)
#remove batch effect
cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "group")


#reduce dimensions
cds <- reduce_dimension(cds)

#determine how many dimensions to use
pdf("00_dim.pdf")
plot_pc_variance_explained(cds)
dev.off()

#check batch
pdf("00_batch.pdf")
plot_cells(cds, color_cells_by="group", label_cell_groups=FALSE)
dev.off()

pdf("01_UMAP.pdf", width = 7, height = 7)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "Idents.BRAIN.",group_label_size = 5,cell_size = 1, label_cell_groups=FALSE)
dev.off()

#clustering
cds <- cluster_cells(cds)
pdf("02_partition.pdf", width = 7, height = 7)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "partition")
dev.off()

#learn graph
cds <- learn_graph(cds)
pdf("03_trajectorise.pdf", width = 7, height = 4)
plot_cells(cds,
           color_cells_by = "Idents.BRAIN.",
           label_groups_by_cluster=TRUE,
           cell_size = 1,
           label_leaves=FALSE,
           label_branch_points=FALSE,
		   label_cell_groups=FALSE)+
           scale_colour_discrete(type=mypalette)
dev.off()

pdf("03_trajectorise_phase.pdf", width = 4, height = 4)
plot_cells(cds,
           color_cells_by = "phase",
           label_groups_by_cluster=TRUE,
           cell_size = 1,
           label_leaves=FALSE,
           label_branch_points=FALSE,
		   label_cell_groups=FALSE)
dev.off()

pdf("03_trajectorise_cell.type.pdf", width = 7, height = 7)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=TRUE,
           cell_size = 0.7,
           label_leaves=FALSE,
           label_branch_points=FALSE,
		   label_cell_groups=FALSE)
dev.off()

get_earliest_principal_node1 <- function(cds, time_bin=opt$node1){

  cell_ids <- which(colData(cds)[,"cell.type"] == time_bin)
 
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]
 
  root_pr_nodes
}

get_earliest_principal_node2 <- function(cds, time_bin=opt$node2){

  cell_ids <- which(colData(cds)[,"cell.type"] == time_bin)

  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}


cds = order_cells(cds, root_pr_nodes=c(get_earliest_principal_node1(cds),get_earliest_principal_node2(cds)))

pdf("trajectorise_set_root.pdf", width = 6, height = 4)
plot_cells(cds,
           color_cells_by = "pseudotime",
		   cell_size = 1,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()

pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
BRAIN@meta.data$pseudotime<-pseudotime


#change the first name of images
b<-paste0("image_",gsub("-",".",names(table(BRAIN$group))))
c<- names(BRAIN@images)
names(BRAIN@images)[1]<-b[!b%in%c]


library(gridExtra)
p<-list()
for(i in names(table(BRAIN$group)))
{
imagename<-paste0("image_",gsub("-",".",i))
temp<-cbind(BRAIN@meta.data[BRAIN@meta.data$group==i,c('group','pseudotime')],BRAIN@images[[imagename]]@coordinates[,1:2])
#h <- as.numeric(max(temp$y) - min(temp$y) + 1)
#w <- as.numeric(max(temp$x) - min(temp$x) + 1)
p[[i]]<-ggplot(temp) + geom_point(mapping=aes(x=x, y=y, colour = pseudotime),size = 6)+
  scale_color_viridis(option="plasma",end=max(temp$pseudotime)/max(BRAIN@meta.data$pseudotime))+
  theme_bw() + theme(panel.grid=element_blank())+
  ggtitle(i)
}
sum<-grid.arrange(grobs=p,ncol=4)
ggsave2("pseudotime_spatial.pdf",sum,width=25,height=20)


#finding genes changing along pseudotime
#The data frame pr_graph_test_res has the Moran's I test results for each gene in the cell_data_set. 
#If you'd like to rank the genes by effect size, sort this table by the morans_Icolumn, which ranges from -1 to +1. 
#A value of 0 indicates no effect, while +1 indicates perfect positive autocorrelation and suggests that nearby cells have very similar values of a gene's expression. 
#Significant values much less than zero are generally rare.
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph=opt$graph, cores=4)
write.csv(file="top_genes.csv",ciliated_cds_pr_test_res)
top_genes<-ciliated_cds_pr_test_res %>% top_n(n = 12, wt = morans_I)
pdf("top_genes.pdf")
plot_cells(cds, genes=top_genes$gene_short_name,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
dev.off()


for (i in top_genes$gene_short_name)
{
library(gridExtra)
p<-list()
for(j in b)
{
p[[j]]<- SpatialFeaturePlot(BRAIN,images=j,features=i,stroke=0,pt.size=12)+ggtitle(j)+scale_color_viridis(option="plasma")
}
sum<-grid.arrange(grobs=p,ncol=4)
ggsave2(paste0(i,"sub_dimspatial.pdf"),sum,width=20,height=20)

}

pr_deg <- subset(ciliated_cds_pr_test_res, q_value < 0.01)
pr_deg_ids<-rownames(pr_deg)
write.csv(file="pr_deg.csv",pr_deg,sep=",",quote=FALSE)

#collect the trajectory-variable genes into modules:
set.seed(1234)
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-3)
write.table(gene_module_df, "gene_module_df.txt",col.names = F)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$Idents.BRAIN.)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

#pdf("pr_deg_heatmap.pdf")
plot1<-pheatmap::pheatmap(agg_mat,scale="row", clustering_method="ward.D",color=plasma(14))
#plot1
#dev.off()
ggsave2("pr_deg_heatmap.pdf",plot1,width=3,height=5)

write.csv(file="agg_mat.csv",agg_mat,sep=",",quote=FALSE)

pr_deg$id<-pr_deg$gene_short_name
temp<-left_join(gene_module_df,pr_deg,by='id')
write.csv(file="final_gene.csv",temp,row.names=F,quote=F)
#plotting important genes



target_gene<-gene_module_df[gene_module_df$module==15,'id']
target_pr<-pr_deg[pr_deg$gene_short_name %in% target_gene$id,]
target_pr[order(target_pr$morans_I,decreasing = T),]%>%head(n=10)

target_genes<-"OL16419,OL28852,OL29365"
target_genes<-unlist(strsplit(target_genes,","))

for (i in target_genes)
{
library(gridExtra)
p<-list()
for(j in b)
{
p[[j]]<- SpatialFeaturePlot(BRAIN,images=j,features=i,stroke=0,pt.size=12,alpha = c(0.1, 1),)+ggtitle(j)
}
sum<-grid.arrange(grobs=p,ncol=4)
ggsave2(paste0(i,"sub_dimspatial.pdf"),sum,width=20,height=20)

}

lineage_cds <- cds[rowData(cds)$gene_short_name %in% target_genes,]
pdf("selected_genes_pseudotime.pdf",width=7,height=4)
plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="Idents.BRAIN.",
                         min_expr=0.01)+
scale_colour_discrete(type=mypalette)
dev.off()

pdf("selected_genes_UMAP.pdf",width=7,height=2.5)
plot_cells(cds, genes=target_genes,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
dev.off()



if(FALSE)
{
cell_wall_genes<-c("OL15002","OL16110","OL19485","OL20442","OL21534","OL24678","OL29397","OL32289")
cell_wall_genes2<-"OL03407,OL06358,OL10532,OL11651,OL15703,OL18437,OL20224,OL23091,OL24482"
cell_wall_genes2<-unlist(strsplit(cell_wall_genes2,","))
cell_wall_genes3<-"OL11924,OL13087,OL13611,OL15091,OL16214,OL17249"
cell_wall_genes3<-unlist(strsplit(cell_wall_genes3,","))
cell_cycle<-"OL00594,OL00947,OL01683,OL04346,OL05207,OL06404,OL07463,OL08028,OL09624,OL10847,OL12726,OL13993,OL15860,OL16920,OL20373,OL22844,OL23098,OL25925,OL27892,OL29689,OL32035"
cell_cycle<-unlist(strsplit(cell_cycle,","))
lineage_cds <- cds[rowData(cds)$gene_short_name %in% cell_wall_genes,]
pdf("cell_wall_genes.pdf",height=14)
plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="phase",
                         min_expr=0.5)
dev.off()

gene_module<-data.frame(id=c(cell_wall_genes,cell_wall_genes3,cell_wall_genes2),
                        module=c(rep("Cell wall organization or biogenesis",length(cell_wall_genes)),
						         rep("Cell wall organization or biogenesis 3",length(cell_wall_genes3)),
								 rep("Cell wall organization or biogenesis 2",length(cell_wall_genes2))))
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group_df)
pdf("selected_go.pdf",height=3)
pheatmap::pheatmap(agg_mat,scale="row", cluster_cols = FALSE,clustering_method="ward.D2")
dev.off()
}
