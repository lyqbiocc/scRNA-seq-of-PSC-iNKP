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
library(corplot)

NK_iNK_integ=readRDS("./WT_tr_iNK_integ_umap.rds")
NK_iNK_integ$group=gsub("_iNK","_CD19CAR-R4iNK",NK_iNK_integ$group)
NK_iNK_integ$group=gsub("_aNK","_activated NK",NK_iNK_integ$group)
NK_iNK_integ$group=gsub("_rNK","_resting NK",NK_iNK_integ$group)
NK_iNK_integ$group[NK_iNK_integ$group=="ESC_CD19CAR-R4iNK"]="PSC_iNK"


NK_iNK_integ$group=factor(NK_iNK_integ$group,levels=c("UCB_activated NK","PB_activated NK","PSC_iNK","UCB_resting NK","PB_resting NK","PB_CD19CAR-R4iNK","BM_CD19CAR-R4iNK","liver_CD19CAR-R4iNK","lung_CD19CAR-R4iNK","spleen_CD19CAR-R4iNK"))

p=DimPlot(NK_iNK_integ,split.by = "group", reduction = "umap", label =F,label.size = 5,group.by = "group1",raster = F,pt.size = 0.0001) &scale_color_manual(values=c(colors[c(11:15,6:10)]))&theme_bw()& theme(legend.position = "bottom")&
  guides(color = guide_legend(nrow = 1, byrow = TRUE, 
                              override.aes = list(size = 2)))&labs(title="")
p
#ggsave(p,filename = "./res_0722/UMAP_split_WTNK_triNK.pdf",width = 11,height = 2.5)

## correlation analysis  
all.genes=rownames(NK_iNK_integ)[rowSums(NK_iNK_integ@assays$RNA@counts>0) > 30]




expr=AverageExpression(NK_iNK_integ,slot = "data",features = all.genes,group.by = "group")
cor_value=round(cor(expr$RNA[,1:5],expr$RNA[,6:10],method = "pearson"),3)
head(cor_value)




library(RColorBrewer)
library(corrplot)
color_1<-rev(brewer.pal(10,"RdYlGn"))
color_2<-colorRampPalette(c("CornflowerBlue","White","tomato"))  #形成一个连续型颜色函数

#pdf(file="./res_0722/cor_NK_heatmap.pdf",width = 4.5,height = 4)
p=corrplot(cor_value,col = color_1,addCoef.col = "black",is.corr = F,col.lim = c(0.5,1),tl.col = "black",cl.align.text = "l",cl.offset = 0.2)
p
#corHybridPlot <- recordPlot()
#dev.off()
