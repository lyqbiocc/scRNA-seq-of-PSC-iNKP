# Figure S8B-S8C: (B) UMAP plots showing induced cells annotated by the Seurat reference mapping workflow with scRNA-seq data set from Rebuffet et al. Each dot represents one cell, and colors represent one iNK cell identity. Plot is faceted by natural NK cell identities. 
# (C) UMAP plots showing induced cells annotated by the Seurat reference mapping workflow with scRNA-seq data set from Jaeger et al. Each dot represents one cell, and colors represent one iNK cell identity. Plot is faceted by natural NK/ILC cell identities.
# library packages  
library(Seurat)

# load data  
## PSC derived-iNKP recipients' iNK cell in vivo  
invivo_iNK <- readRDS("./iNK_integreted.rds")
## natural ILC cells scRAN-seq data was downloaded from GEO databaase (accession number: GSE240441) (PMID: 38956380)   
ILC <- readRDS("./AllTissues_Fig5a_b_c.rds")

invivo_iNK$tissue<-paste0(invivo_iNK$tissue,"_iNK")
ILC$tissue <- paste0(ILC$tissue,"_NK")

# normalization and projection  
ILC <- SCTransform(ILC)
invivo_iNK<- SCTransform(invivo_iNK)
anchorF <- FindTransferAnchors(reference = ILC,query = invivo_iNK, normalization.method = "SCT",reference.reduction = "pca",dims = 1:30)

map_iNK_ILC <- MapQuery(anchorset = anchorF,
         query = invivo_iNK, 
         reference = ILC, 
         refdata = ILC$tissue_chief)

DimPlot(map_iNK_ILC,group.by = "tissue",split.by = "predicted.id",reduction = "umap")
saveRDS(map_iNK_ILC,"./map_ILC_iNK.rds")

library(anndata)
## natural NK cells scRAN-seq data was downloaded from an interactive portal (https://collections.cellatlas.io/meta-nk), the interactive portal was reported by the research (PMID: 38956378)  
tis_nk<- read_h5ad("./PBMC_Noreg.h5ad")
tis_nk_seu <- CreateSeuratObject(counts = t(tis_nk$X), meta.data = tis_nk$obs)
features <- SelectIntegrationFeatures(object.list = list(invivo_iNK,tis_nk_seu), nfeatures = 2000)

#iNK_tis_list <- PrepSCTIntegration(object.list = list(invivo_iNK,tis_nk_seu), anchor.features = features)
anchorF <- FindTransferAnchors(reference = tis_nk_seu,query = invivo_iNK, normalization.method = "LogNormalize",reference.reduction = "pca",dims = 1:30)
map_iNK <- MapQuery(anchorset = anchorF,
                    query = invivo_iNK, 
                    reference = tis_nk_seu, 
                    refdata = tis_nk_seu$ident)

DimPlot(map_iNK,group.by = "tissue",split.by = "predicted.id",reduction = "umap")
saveRDS(map_iNK,"./map_tis_iNK.rds")


