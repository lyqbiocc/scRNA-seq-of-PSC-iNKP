# scRNA-seq of D16-D20 PSC-derived cells
## Figure supplementary 1 A-B
# library packages

library(ggplot2)
#library(ggrepel)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(magrittr)
library(pheatmap)
library(clusterProfiler)
library(cowplot) #combine plots

# D16-D20 PSC-derived cells analysis  
## read 10x data after aggregation

data.seu <- Read10X(data.dir = "./aggr/count/filtered_feature_bc_matrix/")
meta=read.delim("./aggr/aggregation.csv",header = T,sep=",")

data_meta=data.frame(bcnum=colnames(data.seu),reshape2::colsplit(colnames(data.seu), "[-]", names = c("bc","num")),group="")
table(data_meta$num)

data_meta[data_meta$num==1,]$times=meta$sample_id[1]
data_meta[data_meta$num==2,]$times=meta$sample_id[2]
data_meta[data_meta$num==3,]$times=meta$sample_id[3]
data_meta[data_meta$num==4,]$times=meta$sample_id[4]
data_meta[data_meta$num==5,]$times=meta$sample_id[5]

head(data_meta)
rownames(data_meta)=data_meta$bcnum

D16_D20.seu <- CreateSeuratObject(counts = data.seu, project = "GFPneg", min.cells = 5, min.features = 100,meta.data  = data_meta)
D16_D20.seu
# QC  
# global expression of scRNA-seq 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
D16_D20.seu[["percent.mt"]] <- PercentageFeatureSet(D16_D20.seu, pattern = "^MT-")
P=VlnPlot(D16_D20.seu, group.by = "times",features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 4,pt.size = 0,raster = F)&labs(x="",y="")
P

D16_D20.seu<- subset(D16_D20.seu, subset=nFeature_RNA>1000&nFeature_RNA<10000&percent.mt<10)
table(D16_D20.seu$times)

D16_D20.seu=NormalizeData(D16_D20.seu,normalization.method = "LogNormalize", scale.factor =median(D16_D20.seu$nCount_RNA))

## clustering and UMAP  

D16_D20.seu <- FindVariableFeatures(D16_D20.seu, selection.method = "vst", nfeatures = 2000)
D16_D20.seu <- ScaleData(D16_D20.seu, features = VariableFeatures(object = D16_D20.seu))
D16_D20.seu <- RunPCA(D16_D20.seu, features = VariableFeatures(object = D16_D20.seu))
ElbowPlot(D16_D20.seu,ndims = 50)
xx <- cumsum(D16_D20.seu[["pca"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.9)
## do UMAP, and dimmension plot for cell clusters  
# run umap
gc()
## parameter 'dims': the number of dimension of PCA used to umap analysis
D16_D20.seu=FindNeighbors(D16_D20.seu, dims = 1:30)
D16_D20.seu=FindClusters(D16_D20.seu, resolution = 0.3, algorithm = 1)
D16_D20.seu <- RunUMAP(D16_D20.seu, dims = 1:30,return.model=T)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# plot the UMAP result of all cells with DimPlot
## parameter 'pt.size': the point size in UMAP plot
p=plot_grid(
  DimPlot(D16_D20.seu, reduction = "umap", label = T,label.size = 5,group.by = "RNA_snn_res.0.1",raster = F),
  DimPlot(D16_D20.seu, reduction = "umap", label = T,label.size = 5,group.by = "RNA_snn_res.0.2",raster = F),
  DimPlot(D16_D20.seu, reduction = "umap", label = T,label.size = 5,group.by = "RNA_snn_res.0.3",raster = F),
  DimPlot(D16_D20.seu, reduction = "umap", label = T,label.size = 5,group.by = "times",raster = F),
  ncol = 3
)
p & theme_bw()

colors= c('#00f000','#0000f0','#b30000','#f0f000',
          '#00f0f0','#a000f0','#f0a000','#7e3343',
          '#00f0a0','#fb8072','#80b1d3','#fdb462',
          '#b3de69','#fccde5','#bc80bd','#ccebc5',
          '#ffed6f','#64b267','#47a7bd','#f36621',
          '#31629d','#9fde00','#ffbe2a','#ec008c','#ff7404')
## Figure S1B: (B) UMAP visualization of the feature genes expressed in iStroma, iEndo, iHPC, iMkE, iMyeloid, and iNK lineage cells.
p=VlnPlot(D16_D20.seu,features = c("PDGFRA","COL1A2", "CDH5","ECSCR","NR2F2", "PTPRC","CD34","ANGPT1","MEIS1","MECOM","SPINK2", "ITGA2B", "GATA1","GATA2","KLF1","CPA3","CEBPE","CEACAM8","CD68","C1QA","IGLL1", "IL2RB","ITGB7","KIT","IL7R","NKG7","CD38","CD7","MME","IL3RA","IFITM1","KLRC3","KLRD1","TBX21","NCAM1","GNLY"),group.by = "RNA_snn_res.0.3",pt.size = 0,ncol = 8)&scale_fill_manual(values = colors)&theme(axis.text.x = element_text(size = 12))
p

Idents(D16_D20.seu)=D16_D20.seu$RNA_snn_res.0.3
D16_D20.seu=RenameIdents(D16_D20.seu,
                         "7"="iNK lineage",
                         "2"="Hematopoietic_progenitor",
                         "9"="iEndo", 
                         "1"="iStroma",
                         "3"="iMkE",
                         "6"="iMkE",
                         "0"="iMyeloid","4"="iMyeloid","5"="iMyeloid","8"="iMyeloid","10"="iMyeloid")

D16_D20.seu$anno_2=as.vector(Idents(D16_D20.seu))
D16_D20.seu$anno_3=D16_D20.seu$anno_2
saveRDS(D16_D20.seu,file="./subsets_D16_D20_anno3.rds")

# output the main meta information of D16_D20.seu
meta=as.data.frame(D16_D20.seu@meta.data[,c("times","RNA_snn_res.0.3","anno_1","anno_2","anno_3")])
meta$cellID=rownames(D16_D20.seu@meta.data)
write.table(meta,file = "./meta_subsets_D16_D20_anno3.txt",sep="\t",row.names=T,col.names = T)

#  the code for plotting the Day 16 to Day 20 cells
## Figure S1A:(A) UMAP clustering of all the GFPcells of the organoids from days 16, 17, 18, 19, and 20. 
## iStroma, induced Stroma cells. iEndo, induced Endothelial cells. iHPC, induced Hematopoietic progenitor cells. iMkE, induced Megakaryocytic–Erythroid cells, iMyeloid, induced Myeloid cells, iNK lineage cells, induced iNK progenitor cells and iNK cells.
## plot_the_cellpopulation_D16_D20.R






