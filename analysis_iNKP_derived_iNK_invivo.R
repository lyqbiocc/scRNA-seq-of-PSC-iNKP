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

# tissue iNK analysis  
## input aggr 10X data  

data.seu <- Read10X(data.dir = "./aggr/count/filtered_feature_bc_matrix/")
meta=read.delim("./aggr/aggregation.csv",header = T,sep=",")

data_meta=data.frame(bcnum=colnames(data.seu),reshape2::colsplit(colnames(data.seu), "[-]", names = c("bc","num")),group="")
table(data_meta$num)

data_meta[data_meta$num==1,]$group=meta$sample_id[1]
data_meta[data_meta$num==2,]$group=meta$sample_id[2]
data_meta[data_meta$num==3,]$group=meta$sample_id[3]
data_meta[data_meta$num==4,]$group=meta$sample_id[4]
data_meta[data_meta$num==5,]$group=meta$sample_id[5]
#meta[meta$num==19,]$group="D21_GFPneg"

head(data_meta)
rownames(data_meta)=data_meta$bcnum

tr_iNK_aggr <- CreateSeuratObject(counts = data.seu, project = "tissue_iNK", min.cells = 5, min.features = 100,meta.data  = data_meta)
tr_iNK_aggr

rm(data.seu,data_meta)
gc()

# global expression of scRNA-seq 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
tr_iNK_aggr[["percent.mt"]] <- PercentageFeatureSet(tr_iNK_aggr, pattern = "^MT-")

P=VlnPlot(tr_iNK_aggr, group.by = "group",features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 4,pt.size = 0,raster = F)&labs(x="",y="")
P

# global expression of scRNA-seq 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
tr_iNK_aggr<- subset(tr_iNK_aggr, subset=nFeature_RNA<7500&nCount_RNA>500&nCount_RNA<40000&percent.mt<10)
table(tr_iNK_aggr$group)

tr_iNK_aggr=NormalizeData(tr_iNK_aggr,normalization.method = "LogNormalize", scale.factor =median(tr_iNK_aggr$nCount_RNA))

## GC expression  

tr_iNK_aggr$GC_ident=""
tr_iNK_aggr$GC_ident[tr_iNK_aggr@assays$RNA@counts["CXCR4",]>0]="GCpos"
tr_iNK_aggr$GC_ident[tr_iNK_aggr@assays$RNA@counts["CXCR4",]==0]="GCneg"
table(tr_iNK_aggr$group,tr_iNK_aggr$GC_ident)

## subset GC+ cells to analysis  
gc_cells=WhichCells(tr_iNK_aggr,expression = CXCR4>0,slot = "counts")
tr_iNK_aggr_GC=subset(tr_iNK_aggr,cells=gc_cells)
rm(tr_iNK_aggr)
gc()

tr_iNK_aggr_GC$group=gsub("_iNK","_CD19CAR-R4iNK",tr_iNK_aggr_GC$group)

## scale data to do dimension reduction analysis
tr_iNK_aggr_GC <- ScaleData(tr_iNK_aggr_GC, features = VariableFeatures(object = tr_iNK_aggr_GC))
colors= c('#00f000','#0000f0','#b30000','#f0f000',
          '#00f0f0','#a000f0','#f0a000','#7e3343',
          '#00f0a0','#fb8072','#80b1d3','#fdb462',
          '#b3de69','#fccde5','#bc80bd','#ccebc5',
          '#ffed6f','#64b267','#47a7bd','#f36621',
          '#31629d','#9fde00','#ffbe2a','#ec008c','#ff7404')

## gene expression of NK-function genes  

genes=c("NCR1","NCR2","NCR3","NCAM1","FCGR3A","CXCR4", "KLRK1","KLRC1","CD96","KLRD1","FASLG","GZMB","PRF1","TNFSF10","KLRC2","B3GAT1","CD226","CD160","GNLY","IL2RB","KLRC3","KLRF1")

genes=c("NCAM1","FCGR3A","NCR1","NCR3","KLRK1","CD69","KLRC2","KLRC3","SLAMF7","CD244","CD226","KLRB1","KLRC1","KLRD1", "CD96","GZMB","PRF1","GNLY","FASLG","IL2RB")
tr_iNK_aggr_GC$group1=factor(tr_iNK_aggr_GC$group,levels = c( "PB_CD19CAR-R4iNK","BM_CD19CAR-R4iNK", "liver_CD19CAR-R4iNK","lung_CD19CAR-R4iNK","spleen_CD19CAR-R4iNK"))
#colors = c("#E8554E", "#FFD266", "#2AA876",'#a000f0')
p1 = Seurat::VlnPlot(tr_iNK_aggr_GC, features = genes, pt.size = 0,ncol = 10,group.by = "group1")& theme_bw() & theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  strip.text = element_text(face="italic",size = 8,vjust = 0.5), 
  plot.title = element_text(face="italic",size = 15,hjust = 0.5,vjust = 0.5),
  legend.position = "none",
  axis.title.x = element_blank(),
  axis.text.x = element_blank())& labs(y="")&scale_fill_manual(values = colors[c(6:10)])
#scale_fill_manual(values = c("#F8766D","#00B0F6","#00BF7D","#E76BF3"))#&geom_jitter(size=0.001)

p1=p1+theme(legend.position = "bottom",legend.title = element_blank())+guides(size = guide_legend(nrow = 1, byrow = T))+theme(plot.margin = margin(t = 20,  # 顶部边缘距离
                                                                                                                                                   r = 10,  # 右边边缘距离
                                                                                                                                                   b = 10,  # 底部边缘距离
                                                                                                                                                   l = 10)) # 左边边缘距离

p1
#ggsave(p1,filename = "./iNKP_derived_iNK_res/NK_genes_triNK_1012_1.pdf",width = 20,height = 4.2)

## intergration analysis of tr_iNK (invivo) and natural NK cells from PB, UCB and iNK (invitro)  
### input data


# natural NK
PB_rNK=readRDS("../iNK_comp/PB_rNK.rds")
PB_rNK$orig.ident="PB_resting NK"
PB_rNK$group="PB_resting NK"
PB_aNK=readRDS("../iNK_comp/PB_aNK.rds")
PB_aNK$orig.ident="PB_activated NK"
PB_aNK$group="PB_activated NK"
UCB_rNK=readRDS( "../iNK_comp/UCB_rNK.rds")
UCB_rNK$orig.ident="UCB_resting NK"
UCB_rNK$group="UCB_resting NK"
UCB_aNK=readRDS("../iNK_comp/UCB_activated_NK.rds")
UCB_aNK$orig.ident="UCB_activated NK"
UCB_aNK$group="UCB_activated NK"
# iNK
PSC_iNK=readRDS("../iNK_comp/ESC_OL.rds")
PSC_iNK$orig.ident="PSC_iNK"
PSC_iNK$group="PSC_iNK"


NK_iNK_integ=merge(x=tr_iNK_aggr_GC,
                   y=c(NormalizeData(UCB_rNK,normalization.method = "LogNormalize", scale.factor =median(UCB_rNK$nCount_RNA)),
                       NormalizeData(UCB_aNK,normalization.method = "LogNormalize", scale.factor =median(UCB_aNK$nCount_RNA)),
                       NormalizeData(PB_rNK,normalization.method = "LogNormalize", scale.factor =median(PB_rNK$nCount_RNA)),
                       NormalizeData(PB_aNK,normalization.method = "LogNormalize", scale.factor =median(PB_aNK$nCount_RNA)),
                       NormalizeData(PSC_iNK,normalization.method = "LogNormalize", scale.factor =median(PSC_iNK$nCount_RNA)))
)
#rm(tr_iNK_aggr_GC,UCB_rNK,UCB_aNK,PB_rNK,PB_aNK,PSC_iNK)
gc()


### data QC

# global expression of scRNA-seq 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
NK_iNK_integ[["percent.mt"]] <- PercentageFeatureSet(NK_iNK_integ, pattern = "^MT-")
NK_iNK_integ[["percent.ribo"]] <- PercentageFeatureSet(NK_iNK_integ, pattern = "^RP[SL]")
#rb.genes <- rownames(sce)[grep("^RP[SL]",rownames(sce))]
#percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
#sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")

NK_iNK_integ[["log10GenesPerUMI"]] <- log10(NK_iNK_integ$nFeature_RNA)/log10(NK_iNK_integ$nCount_RNA)
NK_iNK_integ[["UMIsPerGene"]] <- NK_iNK_integ$nCount_RNA/NK_iNK_integ$nFeature_RNA

# erythrocyte-related genes-hemoglobin family
##mouse: HB.genes <- "Hba-a1/Hba-a2/Hbb/Hbd/Hbb-y/Hbb-bh1/hbb-bh2/Hbb-b1/Hbb-b2/hbb-ar/hbb-bt/hbb-bs/lrp5/Hbq1a/Hbq1b/hba/Hba-x" %>% strsplit(.,split = "/") %>% unlist %>% capitalize(.)
HB.genes <- "HBA1/HBA2/HBB/HBD/HBE1/HBG1/HBG2/HBM/HBQ1/HBZ" %>% strsplit(.,split = "/") %>% unlist %>% toupper(.) %>% intersect(.,rownames(NK_iNK_integ))
if(!is.null(HB.genes)){
  percent.hb <- Matrix::colSums(NK_iNK_integ@assays$RNA@counts[HB.genes,])/Matrix::colSums(NK_iNK_integ@assays$RNA@counts)*100
  NK_iNK_integ <- AddMetaData(NK_iNK_integ, percent.hb, col.name = "percent.hb")
}

# global expression of scRNA-seq 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#NK_iNK_integ[["percent.mt"]] <- PercentageFeatureSet(NK_iNK_integ, pattern = "^MT-")

#NK_iNK_integ$group=factor(NK_iNK_integ$group,levels = c("PB_rNK","PB_aNK","UCB_rNK","UCB_aNK","ESC_iNK","PB_iNK","BM_iNK","liver_iNK","lung_iNK","spleen_iNK"))
P=VlnPlot(NK_iNK_integ, group.by = "group",features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.ribo","percent.hb","log10GenesPerUMI","UMIsPerGene"), ncol = 4,pt.size = 0,raster = F)&labs(x="",y="")&theme_bw()&NoLegend()&scale_fill_manual(values = colors[c(11:15,6:10)])&theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),plot.title   = element_text(hjust=0.5,vjus=1))
P



NK_iNK_integ=subset(NK_iNK_integ,subset=percent.mt<15)


### clustering 
NK_iNK_integ <- FindVariableFeatures(NK_iNK_integ, selection.method = "vst", nfeatures = 2000)
NK_iNK_integ <- ScaleData(NK_iNK_integ, features = VariableFeatures(object = NK_iNK_integ))
NK_iNK_integ=RunPCA(NK_iNK_integ, features = VariableFeatures(object = NK_iNK_integ))


ElbowPlot(NK_iNK_integ,ndims = 50)
xx <- cumsum(NK_iNK_integ[["pca"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.9)


### cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# CellCycleScoring() parameter
## parameter 's.features': a vector of genes associated with S phase
## parameter 'g2m.features': a vector of genes associated with G2M phase
## parameter 'set.ident': if true, sets identity to phase assignments

# REMOVE all cell cycle effect
## G2M.Score & S.Score<0 则定义为G1期
NK_iNK_integ <- CellCycleScoring(NK_iNK_integ, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

NK_iNK_integ$CC.Difference <- NK_iNK_integ$S.Score - NK_iNK_integ$G2M.Score
#GFPneg_aggr <- ScaleData(GFPneg_aggr, vars.to.regress = "CC.Difference", features = rownames(GFPneg_aggr))
NK_iNK_integ <- ScaleData(NK_iNK_integ, vars.to.regress = c("group","S.Score","G2M.Score","CC.Difference"), features = VariableFeatures(object = NK_iNK_integ))


### do UMAP, and dimmension plot for cell clusters  

# run umap
gc()
## parameter 'dims': the number of dimension of PCA used to umap analysis
NK_iNK_integ=FindNeighbors(NK_iNK_integ, dims = 1:40)#1:40
NK_iNK_integ=FindClusters(NK_iNK_integ, resolution = 0.1, algorithm = 1)
NK_iNK_integ=FindClusters(NK_iNK_integ, resolution = 0.2, algorithm = 1)
NK_iNK_integ=FindClusters(NK_iNK_integ, resolution = 0.3, algorithm = 1)
NK_iNK_integ=FindClusters(NK_iNK_integ, resolution = 0.4, algorithm = 1)
NK_iNK_integ <- RunUMAP(NK_iNK_integ, dims = 1:40,return.model=T)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# plot the UMAP result of all cells with DimPlot
## parameter 'pt.size': the point size in UMAP plot
p=plot_grid(
  DimPlot(NK_iNK_integ, reduction = "umap", label = T,label.size = 5,group.by = "RNA_snn_res.0.1",raster = F),
  DimPlot(NK_iNK_integ, reduction = "umap", label = T,label.size = 5,group.by = "RNA_snn_res.0.2",raster = F),
  DimPlot(NK_iNK_integ, reduction = "umap", label = T,label.size = 5,group.by = "RNA_snn_res.0.3",raster = F),
  DimPlot(NK_iNK_integ, reduction = "umap", label = T,label.size = 5,group.by = "RNA_snn_res.0.4",raster = F),
  DimPlot(NK_iNK_integ, reduction = "umap", label = T,label.size = 5,group.by = "group",raster = F),
  DimPlot(NK_iNK_integ, reduction = "umap", label = T,label.size = 5,group.by = "Phase",raster = F),
  ncol = 3
)
p & theme_bw()

saveRDS(NK_iNK_integ,file = "./WT_tr_iNK_integ_umap.rds")
