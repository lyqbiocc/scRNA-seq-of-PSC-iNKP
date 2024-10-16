# D16-D20 PSC-derived cells analysis  
## plot the cell population from D16-D20 PSC-derived cells
#D16_D20.seu=readRDS("./subsets_D16_D20_anno3.rds")
D16_D20.seu$anno_2_f=factor(D16_D20.seu$anno_2,levels = c("iStroma","iEndo","Hematopoietic_progenitor","iMkE","iMyeloid","iNK lineage"))
p=DimPlot(D16_D20.seu,pt.size = 0.001,raster = F,group.by = c("anno_2_f"),label = F,label.size = 0,ncol = 1)+scale_fill_brewer(palette = "Dark2")+labs(title="")&theme(
  axis.text.x=element_text(colour="black",size=8),
  #strip.text = element_text(size = 15,color = "black"),
  axis.text.y=element_text(colour="black",size=8),
  axis.title.x=element_text(size = 8),
  axis.title.y=element_text(size = 8),
  panel.border = element_rect(size=1,color = "black"),
  axis.line = element_line(colour = "black",size=0),
  legend.text=element_text( colour="black",size=8),
  legend.title=element_text(colour="black",size=8),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # plot.title = element_rect(fill=NA,size = 1),
  #plot.title = element_text(hjust = 0.5,size=15),
  plot.title = element_blank(),
  legend.position = "bottom")&
  guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 4)))
p
#ggsave(p,filename = "./D16_D20_res/UMAP_lanscap.pdf",width = 5,height = 5)

## plot the lineage-gene expression from D16-D20 PSC-derived cells
p=FeaturePlot(D16_D20.seu,features = c("PTPRC","PDGFRA","CDH5","ANGPT1","SPINK2","IGLL1","GATA1","PF4","GATA2","CPA3","SPI1","MPO","CLEC4A","CD34","CD38","CD7","NKG7","IL2RB","GNLY"),ncol = 7,pt.size = 0.01,cols = c("lightgrey","red"))&theme_bw() &
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
p
#ggsave(p,filename = "./D16_D20_res/UMAP_exp.pdf",width = 21,height = 9)