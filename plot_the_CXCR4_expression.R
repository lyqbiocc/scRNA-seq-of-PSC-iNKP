# expression of CXCR4  

## expression of CXCR4 of iNKP and iNK (invitro)

### load rds 
NC_NK_iNKP_iNK=readRDS("./NC_NK_iNKP_iNK.rds")
#NC_NK: PB_NK;UCB_NK
#iNKP: iNKP of Days-16 to Days20 (invitro)
#iNK: iNK of Days27 (invitro)

data4boxplot=data.frame(Sample = colnames(NC_NK_iNKP_iNK$RNA), 
                        stage=NC_NK_iNKP_iNK$group1,
                        CXCR4=GetAssayData(NC_NK_iNKP_iNK,slot = "data")["CXCR4",]
                        
) %>% reshape2::melt()

data4boxplot$stage=factor(data4boxplot$stage,levels = c("PB_NK","UCB_NK","iNKP","iNK"))
box_p2=ggplot(data4boxplot, aes(x = variable, y = value)) +   
  geom_boxplot(aes(fill=stage),outlier.size = 0,outlier.shape = NA)+
  #geom_jitter(width = 0.4, size = 0, alpha = 0.4) + 
  #geom_boxplot(aes(color = stage), width = 0.7, outlier.shape = NA) +
  #stat_boxplot(geom = "errorbar", color = stage, width = 0.3, size = 0.8) +
  scale_colour_manual(values = stage) +
  facet_wrap(.~variable,  scales="free",ncol = 5) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face="italic",size = 8) , legend.position = "bottom",axis.title.x = element_blank(),axis.text.x = element_blank())+labs(x="",y="")#&scale_fill_manual(values = colors[c(15:18)])
#box_p=box_p+theme(legend.position = "bottom",legend.title = element_blank())
box_p2


#BM_NK of GSE149938
NK_cells=readRDS("./NK_cells_GSE149938.rds")
### plot CXCR4 expression  
colors= c('#00f000','#0000f0','#b30000','#f0f000',
          '#00f0f0','#a000f0','#f0a000','#7e3343',
          '#00f0a0','#fb8072','#80b1d3','#fdb462',
          '#b3de69','#fccde5','#bc80bd','#ccebc5',
          '#ffed6f','#64b267','#47a7bd','#f36621',
          '#31629d','#9fde00','#ffbe2a','#ec008c','#ff7404')

data4boxplot=data.frame(Sample = colnames(NK_cells$RNA), 
                        stage=NK_cells$group1,
                        CXCR4=GetAssayData(NK_cells,slot = "data")["CXCR4",]
                        
) %>% reshape2::melt()


#data4boxplot=data4boxplot[data4boxplot$stage=="NKP"|data4boxplot$stage=="NK",]
box_p1=ggplot(data4boxplot, aes(x = variable, y = value)) +   
  geom_boxplot(aes(fill=stage),outlier.size = 0,outlier.shape = NA)+
  #geom_jitter(width = 0.4, size = 0, alpha = 0.4) + 
  #geom_boxplot(aes(color = stage), width = 0.7, outlier.shape = NA) +
  #stat_boxplot(geom = "errorbar", color = stage, width = 0.3, size = 0.8) +
  scale_colour_manual(values = stage) +
  facet_wrap(.~variable,  scales="free",ncol = 5) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face="italic",size = 8) , legend.position = "bottom",axis.title.x = element_blank(),axis.text.x = element_blank())+labs(x="",y="")+scale_fill_manual(values =colors[20:24])
#box_p=box_p+theme(legend.position = "bottom",legend.title = element_blank())
box_p1