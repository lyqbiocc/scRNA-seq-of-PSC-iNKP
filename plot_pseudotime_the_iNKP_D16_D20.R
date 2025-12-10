# pseudo-time plot  
## Figure 1C:  (C) pseudo-time analysis of other-iHPC, e-iNKP, l-iNKP, and im-iNK. 
CD45_CD34_imNK.cds=readRDS("./pesud_CD45_CD34_imNK_1013.rds")
p=plot_cells(CD45_CD34_imNK.cds,cell_size = 0.5,
             color_cells_by = "anno_6_f",
             label_cell_groups=FALSE,
             label_leaves=TRUE,
             label_branch_points=T,
             graph_label_size=1.5,trajectory_graph_segment_size =0.5)+
  scale_color_manual(values =c( "lightgrey",RColorBrewer::brewer.pal(9,"Set1")[c(5,4,8)]))+
  theme_bw()+
  theme(plot.title = element_text(vjust=0,hjust=0.5),
        legend.title = element_text(size=0),
        legend.text = element_text(size=6),
        axis.title.x  = element_text(size=1),
        axis.title.y  = element_text(size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(colour=guide_legend(ncol = 1,bycol= T,override.aes = list(size=1),order = 1))+
  theme(legend.spacing.x = unit(0.1, "cm"),  # 设置水平间距为0.5cm
        legend.spacing.y = unit(0.1, "cm"))#+scale_x_continuous(limits = c(-8,8),breaks = c(-8,-4,0,4,8))+scale_y_continuous(limits = c(-5,6),breaks=c(-5,-2.5,0,3,6))
p
#ggsave(p,filename = "./D16_D20_res/pesudo-time_iNKP.pdf",width = 4.5,height = 3)

