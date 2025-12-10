# Figure 1D: (D) Violin plots of gene expression levels of CD34, CD38, CD7, CD123-encoding IL3RA, CD127-encoding IL7R, CD10-encoding MME, CD56-encoding NCAM1, and CD16-encoding FCGR3A in other-iHPC, iNKP (e-iNKP and l-iNKP) and im-iNK. 
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

CD45_CD34_imNK=readRDS(file="./CD45_CD34_imNK_1016.rds")
#CD45_CD34_imNK$anno_3=as.vector(Idents(CD45_CD34_imNK))
CD45_CD34_imNK$anno_7=as.vector(Idents(CD45_CD34_imNK))
CD45_CD34_imNK$anno_7[CD45_CD34_imNK$anno_7=="e-iNKP"]="iNKP"
CD45_CD34_imNK$anno_7[CD45_CD34_imNK$anno_7=="l-iNKP"]="iNKP"
CD45_CD34_imNK$anno_7_f=factor(CD45_CD34_imNK$anno_7,levels = c("other_iHPC","iNKP","im-iNK"))

p=VlnPlot(CD45_CD34_imNK,features = c("CD34","CD38", "CD7","IL3RA","IL7R", "MME","NCAM1","FCGR3A"),
          group.by = "anno_7_f",pt.size = 0.001,ncol = 8) & 
scale_fill_manual(values = colors)&theme(axis.text.x = element_text(size = 12))
