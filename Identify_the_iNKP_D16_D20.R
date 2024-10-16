#subset Hematopoietic_progenitor and iNK lineage cells to second-analysis  

D16_D20.seu <- readRDS("subsets_D16_D20_anno3.rds")

CD45_CD34_imNK=subset(D16_D20.seu,cells=c(
  WhichCells(subset(D16_D20.seu,subset=anno_2!="iStroma" &anno_2!="iEndo"),expression =PTPRC>1&CD34>1,slot = "counts"),
  WhichCells(subset(D16_D20.seu,subset=anno_2=="iNK lineage"))))

## identify the highly various genes


all.features=rownames(CD45_CD34_imNK)[(rowSums(CD45_CD34_imNK@assays$RNA@counts>0)>6)] #cell_num*0.05%
CD45_CD34_imNK=subset(CD45_CD34_imNK,features=all.features)
CD45_CD34_imNK <- FindVariableFeatures(CD45_CD34_imNK, selection.method = "vst", nfeatures = 2000)



## scale data to do dimension reduction analysis

CD45_CD34_imNK <- ScaleData(CD45_CD34_imNK, features = VariableFeatures(object = CD45_CD34_imNK))

## do PCA analysis
#backgroud gene list is the HVGs, which can obtain by VariableFeatures()  

CD45_CD34_imNK <- RunPCA(CD45_CD34_imNK, features = VariableFeatures(object = CD45_CD34_imNK))

ElbowPlot(CD45_CD34_imNK,ndims = 50)
xx <- cumsum(CD45_CD34_imNK[["pca"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.9)


## do UMAP, and dimmension plot for cell clusters  

## parameter 'dims': the number of dimension of PCA used to umap analysis
CD45_CD34_imNK=FindNeighbors(CD45_CD34_imNK, dims = 1:30)
CD45_CD34_imNK=FindClusters(CD45_CD34_imNK, resolution = 0.1, algorithm = 1)
CD45_CD34_imNK=FindClusters(CD45_CD34_imNK, resolution = 0.2, algorithm = 1)
CD45_CD34_imNK=FindClusters(CD45_CD34_imNK, resolution = 0.3, algorithm = 1)
CD45_CD34_imNK=FindClusters(CD45_CD34_imNK, resolution = 0.4, algorithm = 1)
CD45_CD34_imNK=FindClusters(CD45_CD34_imNK, resolution = 0.5, algorithm = 1)
CD45_CD34_imNK=FindClusters(CD45_CD34_imNK, resolution = 0.6, algorithm = 1)
CD45_CD34_imNK <- RunUMAP(CD45_CD34_imNK, dims = 1:30,return.model=T)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# plot the UMAP result of all cells with DimPlot
## parameter 'pt.size': the point size in UMAP plot
p=plot_grid(
  DimPlot(CD45_CD34_imNK, reduction = "umap", label = T,label.size = 5,group.by = "RNA_snn_res.0.1",raster = F),
  DimPlot(CD45_CD34_imNK, reduction = "umap", label = T,label.size = 5,group.by = "RNA_snn_res.0.6",raster = F),
  DimPlot(CD45_CD34_imNK, reduction = "umap", label = T,label.size = 5,group.by = "RNA_snn_res.0.5",raster = F),
  DimPlot(CD45_CD34_imNK, reduction = "umap", label = T,label.size = 5,group.by = "RNA_snn_res.0.4",raster = F),
  DimPlot(CD45_CD34_imNK, reduction = "umap", label = T,label.size = 5,group.by = "anno_2",raster = F),
  DimPlot(CD45_CD34_imNK, reduction = "umap", label = T,label.size = 5,group.by = "times",raster = F),
  ncol = 3
)
p & theme_bw()

### locate the iNKP-like cells  
p=FeaturePlot(CD45_CD34_imNK,features = c("PTPRC","CD34","CD38","CD7","MME","IL2RB","KLRC3","NKG7","CD247","IGLL1","CXCR4"))
p

p=DimPlot(CD45_CD34_imNK, reduction = "umap", label = T,label.size = 3,group.by = "anno_2",raster = F,split.by = "times",ncol = 5)
p

Idents(CD45_CD34_imNK)=CD45_CD34_imNK$RNA_snn_res.0.6
CD45_CD34_imNK=RenameIdents(CD45_CD34_imNK,
                            "4"="e-iNKP",
                            "6"="e-iNKP",
                            "2"="l-iNKP", 
                            "0"="im-iNK",
                            "1"="other-iHPC",
                            "3"="other-iHPC",
                            "5"="other-iHPC","7"="other-iHPC","8"="other-iHPC","9"="other-iHPC","10"="other-iHPC","11"="other-iHPC")

CD45_CD34_imNK$anno_3=as.vector(Idents(CD45_CD34_imNK))


## mapping with natural NKP and NK cells  

### load natural cells gene expression (STRT-seq)

natural_immcell=read.delim("./GSE149938_umi_matrix.csv",sep=",",header = T)

dim(natural_immcell)

#rownames(natural_immcell)[grep(rownames(natural_immcell),pattern = "NKP")]
meta=reshape2::colsplit(rownames(natural_immcell),"[_]",names=c("group","orig","bc"))
NC_cells <- CreateSeuratObject(counts = t(natural_immcell),meta.data = meta,project = "natural_blood_cells", min.cells = 5, min.features = 100)
#rm(data.seu)
gc()

# global expression of scRNA-seq 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
NC_cells[["percent.mt"]] <- PercentageFeatureSet(NC_cells, pattern = "^MT-")

P=VlnPlot(NC_cells, group.by = "group",features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0,raster = F)&labs(x="",y="")

P

Idents(NC_cells)=NC_cells$group
NK_cells=subset(NC_cells,idents=c("NKP","kineNK","toxiNK"))
rm(natural_immcell,meta,NC_cells)
gc()

NK_cells=NormalizeData(NK_cells,normalization.method = "LogNormalize", scale.factor =median(NK_cells$nCount_RNA))


### mapping with natural NKP and NK cells  

#### sctransform before IntegrateData  (finalized fig)  
# run SCTransform on each object separately
iNKP_obj=SCTransform(CD45_CD34_imNK, verbose = FALSE)

features1=rownames(NK_cells)[rowSums(NK_cells@assays$RNA@counts>0)>0]
NK_cells=subset(NK_cells,features=features1)
#NK_cells=subset(NK_cells,subset=group1=="NKP"|group1=="NK")
NK_cells=SCTransform(NK_cells, verbose = FALSE)
gc()
# set the maximum allowed size of seurat object list
#options(future.globals.maxSize = 5000 * 1024^2)
# select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated
pb.features <- SelectIntegrationFeatures(object.list = list(iNKP_obj,NK_cells), nfeatures = 2000)
pb.list <- PrepSCTIntegration(object.list = list(iNKP_obj,NK_cells), anchor.features = pb.features,verbose = T)
gc()
# identify anchors and integrate the datasets. Commands are identical to the standard workflow, but make sure to set normalization.method = 'SCT'
pb.anchors <- FindIntegrationAnchors(object.list = pb.list, normalization.method = "SCT",anchor.features = pb.features, verbose = FALSE,k.filter = 100)
pb.integrated <- IntegrateData(anchorset = pb.anchors, normalization.method = "SCT", verbose = FALSE)
gc()

DefaultAssay(pb.integrated) <- "integrated"
pb.integrated <- RunPCA(pb.integrated, verbose = FALSE)
ElbowPlot(pb.integrated,ndims = 30)

pb.integrated <- RunUMAP(pb.integrated, dims = 1:25,n.neighbors = 30)


