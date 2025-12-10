# Figure S1C: Violin plots of gene expression levels of CD56-encoding NCAM1, CD16-encoding FCGR3A, GZMA, GZMB, perforin-encoding PRF1, TRAIL-encoding TNFSF10, FAS, and CD178-encoding FASLG in UCB_activated NK, PSC-iNK, and iNKP.
## library
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

# load data  
## UCB-NK (primary)
UCB_aNK=readRDS("../iNK_comp/UCB_activated_NK.rds")
UCB_aNK$orig.ident="UCB_activated NK"
UCB_aNK$group="UCB_activated NK"
## iNK
PSC_iNK=readRDS("../iNK_comp/ESC_OL.rds")
PSC_iNK$orig.ident="PSC_iNK"
PSC_iNK$group="PSC_iNK"
## iNKP
CD45_CD34_imNK=readRDS(file="./CD45_CD34_imNK_1016.rds")
#CD45_CD34_imNK$anno_7=as.vector(Idents(CD45_CD34_imNK))
#CD45_CD34_imNK$anno_7[CD45_CD34_imNK$anno_7=="e-iNKP"]="iNKP"
#CD45_CD34_imNK$anno_7[CD45_CD34_imNK$anno_7=="l-iNKP"]="iNKP"
#CD45_CD34_imNK$anno_7_f=factor(CD45_CD34_imNK$anno_7,levels = c("other_iHPC","iNKP","im-iNK"))
CD45_CD34_imNK$group=CD45_CD34_imNK$anno_7
iNKP=subset(CD45_CD34_imNK,subset=group=="iNKP")
rm(CD45_CD34_imNK)
gc()

##merge data  
integ=merge(x=UCB_aNK,
                   y=c(NormalizeData(PSC_iNK,normalization.method = "LogNormalize", scale.factor =median(PSC_iNK$nCount_RNA)),
                       NormalizeData(iNKP,normalization.method = "LogNormalize", scale.factor =median(iNKP$nCount_RNA))
)




