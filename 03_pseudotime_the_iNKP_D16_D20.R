# pesudo-times analysis  
## library packages
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(ggplot2)
library(magrittr)
library(pheatmap)
library(clusterProfiler)
library(cowplot) #combine plots
library("monocle3",lib.loc="/public/home/lin_yunqing/.conda/envs/RNAseq/lib/R/library")

## obtain data to create cds object [expression: lognormalization of UMIs]  
CD45_CD34_imNK=readrds("./CD45_CD34_imNK_1016.rds")
CD45_CD34_imNK.meta <- CD45_CD34_imNK@meta.data
CD45_CD34_imNK.gene <- data.frame(gene_short_name = rownames(CD45_CD34_imNK), row.names = rownames(CD45_CD34_imNK))
CD45_CD34_imNK.counts <- GetAssayData(CD45_CD34_imNK, slot = "data") %>% as.matrix #matrix is required

## create cds object  
CD45_CD34_imNK.cds <- new_cell_data_set(CD45_CD34_imNK.counts,
                                     cell_metadata = CD45_CD34_imNK.meta,
                                     gene_metadata = CD45_CD34_imNK.gene)

rm(CD45_CD34_imNK.meta,CD45_CD34_imNK.gene,CD45_CD34_imNK.counts)
gc()
## normalization/scale/pca/umap in monocle  
## use "data"[lognormalization of UMIs] of seurat object, and not re-normalize and sclale data
CD45_CD34_imNK.cds <- preprocess_cds(CD45_CD34_imNK.cds,
                                  #residual_model_formula_str = "~ nFeature_RNA + nCount_RNA",
                                  #use_genes = VariableFeatures(CD45_CD34_imNK),
                                  method = c("PCA"),
                                  #norm_method="log",
                                  norm_method="none",
                                  #pseudo_count=1,
                                  scaling =F,
                                  num_dim = 10)

CD45_CD34_imNK.cds <- reduce_dimension(CD45_CD34_imNK.cds,preprocess_method = "PCA",
                                    reduction_method = c("UMAP"),
                                    umap.n_neighbors = 15,
                                    umap.min_dist = 0.1,
                                    max_components = 2)

## load UMAP result into CDS object  
cds.embed <- CD45_CD34_imNK.cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pb.integrated, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
CD45_CD34_imNK.cds@int_colData$reducedDims$UMAP <- int.embed

## clusteringr in monocle  
CD45_CD34_imNK.cds <- cluster_cells(CD45_CD34_imNK.cds,resolution=0.0001,k=50,cluster_method = c("louvain"))#resolution=0.0001
unique(partitions(CD45_CD34_imNK.cds))
#plot_cells(CD45_CD34_imNK.cds, color_cells_by = "partition",group_cells_by="partition")
plot_cells(CD45_CD34_imNK.cds, color_cells_by = "anno_6",group_label_size = 3)
plot_cells(CD45_CD34_imNK.cds, color_cells_by = "cluster",group_label_size = 3)
table(colData(CD45_CD34_imNK.cds)$anno_6)
table(CD45_CD34_imNK.cds@clusters$UMAP$clusters)

## search pseudotime trajectory
CD45_CD34_imNK.cds <- learn_graph(CD45_CD34_imNK.cds, 
                               use_partition = F,
                               close_loop = F,
                               learn_graph_control=list(minimal_branch_len=7.5,orthogonal_proj_tip=T,prune_graph=T))
plot_cells(CD45_CD34_imNK.cds,
           color_cells_by = "anno_6",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)



## oder cells to search root => order_cells()    
p=FeaturePlot(pb.integrated,features = c("PTPRC","AVP","ANGPT1","SPINK2","IGLL1","GATA1","GATA2","SPI1","CD34","MECON","HLF","HOXA9","HOXA10"),ncol = 7,pt.size = 0.01,cols = c("lightgrey","red"))&theme_bw() &
  theme(
    axis.text.x=element_text(colour="black",size=8),
    axis.text.y=element_text(colour="black",size=8),
    axis.title.x=element_text(size = 8),
    axis.title.y=element_text(size = 8),
    panel.border = element_rect(),
    axis.line = element_line(colour = "black",size=0),
    legend.text=element_text(colour="black",size=8),
    legend.title=element_text(colour="black",size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face="italic",hjust = 0.5,size=10),legend.position = "bottom")&
  #scale_color_manual(values = c(colors[c(1:15)]))&
  guides(color = guide_legend(nrow = 1, byrow = T, override.aes = list(size = 1)))

### manully selection of root via HSC-related gene expression
CD45_CD34_imNK.cds <- order_cells(CD45_CD34_imNK.cds)
CD45_CD34_imNK.cds@principal_graph_aux@listData$UMAP$root_pr_nodes

saveRDS(CD45_CD34_imNK.cds,file = "./D16_D20_res/pesud_CD45_CD34_imNK_1013.rds")

