# visualization, clustering
pb.integrated=readRDS("./iNKP_NKP_mapping.rds")
plots <- DimPlot(pb.integrated, split.by = c("group"),label = T)
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                     override.aes = list(size = 3)))


#early-NKP(e-iNKP)
#late-nkp(l-iNKP)
#immature-iNK(im-iNK)
pb.integrated$group1[rownames(NK_cells@meta.data)]=NK_cells$group1
#NK_cells$group1[NK_cells$group1=="kineNK"|NK_cells$group1=="toxiNK"]="NK"
pb.integrated$group1[rownames(CD45_CD34_imNK@meta.data)]=CD45_CD34_imNK$anno_3
pb.integrated$group1=factor(pb.integrated$group1,levels = c("other-iHPC","NKP", "NK","e-iNKP","l-iNKP","im-iNK"))

## Figure 1B: (B) UMAP illustrating the projection of NKP, NK with e-iNKP, l-iNKP, and im-iNK. NKP, natural NK progenitor cells, NK, natural NK cells, other-iHPC, other induced hematopoietic cells. e-iNKP, early iNK progenitor cells. l-iNKP, late iNK progenitor cells. im-iNK, immature iNK cells. 
p=DimPlot(pb.integrated, group.by  = c("group1"),label = F)&scale_color_manual(values =c( "lightgrey",RColorBrewer::brewer.pal(9,"Set1")[c(3,2,5,4,8)]))&labs(title="")&theme_bw() &
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.text.y=element_text(colour="black",size=8,face="bold"),
    axis.title.x=element_text(size = 8,face="bold"),
    axis.title.y=element_text(size = 8,face="bold"),
    panel.border = element_rect(),
    axis.line = element_line(colour = "black",size=0),
    legend.text=element_text(face="bold", colour="black",size=8),
    legend.title=element_text(face="bold",colour="black",size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face="bold.italic",hjust = 0.5,size=10),legend.position = "bottom")&
  #scale_color_manual(values = c(colors[c(1:15)]))&
  guides(color = guide_legend(nrow = 1, byrow = T, override.aes = list(size = 3)))
p
#ggsave(p,filename = "./D16_D20_res/iNKP_nkp_UMAPmerg.pdf",width = 3,height = 3)



p=DimPlot(pb.integrated, group.by  = c("group1"),label = F,split.by = "times",ncol = 6)&scale_color_manual(values =c( "lightgrey",RColorBrewer::brewer.pal(9,"Set1")[c(3,2,5,4,8)]))&labs(title="")&theme_bw() &
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.text.y=element_text(colour="black",size=8,face="bold"),
    axis.title.x=element_text(size = 8,face="bold"),
    axis.title.y=element_text(size = 8,face="bold"),
    panel.border = element_rect(),
    axis.line = element_line(colour = "black",size=0),
    legend.text=element_text(face="bold", colour="black",size=8),
    legend.title=element_text(face="bold",colour="black",size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face="bold.italic",hjust = 0.5,size=10),legend.position = "bottom")&
  #scale_color_manual(values = c(colors[c(1:15)]))&
  guides(color = guide_legend(nrow = 1, byrow = T, override.aes = list(size = 3)))
p
#ggsave(p,filename = "./D16_D20_res/iNKP_nkp_UMAP.pdf",width = 12,height = 3)


