library(limma)
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(egg)
library(grid)
library(stringr)
library(ggrepel)
library(ggpp)
library(stringr)
library(scales)
library(cowplot)
library(ggpubr)
library(volcano3D)
library(ggfortify)
library(ggdendro)
library(dendextend)
library(scico)
library(standR)
library(SpatialExperiment)
library(pals)
library(ggpattern)
library(tidyverse)
library(patchwork)
library(org.Hs.eg.db)
library(viridis)
library(clusterProfiler)

## enrichGO 

ego_fibroblast_v_fibrosis_up <-enrichGO(rownames(topT.fibroblast_v_fibrosis[topT.fibroblast_v_fibrosis$logFC > 1 & topT.fibroblast_v_fibrosis$adj.P.Val < 0.05,]), OrgDb=org.Hs.eg.db, ont="BP", keyType = "SYMBOL")
saveRDS(ego_fibroblast_v_fibrosis_up, "IPF_ego_fibroblast_v_fibrosis_up.rds")


## import enrichGO outputs 

IPF_ego_fibroblast_v_fibrosis_up <- 
readRDS("./../within_IPF/IPF_ego_fibroblast_v_fibrosis_up.rds")
IPF_ego_fibrosis_v_neutral_up <- 
  readRDS("./../within_IPF/IPF_ego_fibrosis_v_neutral_up.rds")
IPF_ego_fibroblast_v_neutral_up <- 
  readRDS("./../within_IPF/IPF_ego_fibroblast_v_neutral_up.rds")


NSIP_ego_central_v_inflammatory_up <- 
  readRDS("./../within_NSIP/NSIP_ego_central_v_inflammatory_up.rds")
NSIP_ego_peripheral_v_inflammatory_up <- 
  readRDS("./../within_NSIP/NSIP_ego_peripheral_v_inflammatory_up.rds")
# NSIP_ego_peripheral_v_central_up <- 
#   readRDS("./../within_NSIP/NSIP_ego_peripheral_v_central_up.rds")



CHP_ego_fibrosis_v_neutral_up <- 
  readRDS("./../within_CHP/CHP_ego_fibrosis_v_neutral_up.rds")
CHP_ego_granuloma_v_neutral_up <- 
  readRDS("./../within_CHP/CHP_ego_granuloma_v_neutral_up.rds")
CHP_ego_granuloma_v_inflammatory_up <- 
  readRDS("./../within_CHP/CHP_ego_granuloma_v_inflammatory_up.rds")

## plot merged results
up <- merge_result(list(
  "F.foci vs. Fibrosis" = clusterProfiler::simplify(IPF_ego_fibroblast_v_fibrosis_up, cutoff=0.3),
  "Fibrosis vs. Uninvolved" = clusterProfiler::simplify(IPF_ego_fibrosis_v_neutral_up, cutoff=0.3),
  "F.foci vs. Uninvolved" = clusterProfiler::simplify(IPF_ego_fibroblast_v_neutral_up, cutoff=0.3),
  
  "Central F. vs. Inflam." = clusterProfiler::simplify(NSIP_ego_central_v_inflammatory_up, cutoff=0.3),
  "Peripheral F. vs. Inflam." = clusterProfiler::simplify(NSIP_ego_peripheral_v_inflammatory_up, cutoff=0.3)
  #NSIP_peri_central = NSIP_ego_peripheral_v_central_up,
  
  #CHP_fibrosis_neutral = CHP_ego_fibrosis_v_neutral_up,
  #CHP_granuloma_neutral = CHP_ego_granuloma_v_neutral_up,
  #CHP_granuloma_inflam = CHP_ego_granuloma_v_inflammatory_up
)) %>% dotplot(., showCategory=3)+
  scale_colour_viridis(limits=c(as.double(1e-20), 5e-2),
                       trans="log",
                       breaks=c(1e-20, 1e-10, 1e-5, 5e-2))+
  # scale_size_continuous(limits=c(1,15),
  #                       breaks=c(5,10,15))+
  # 
  # scale_x_continuous(
  #   labels=scales::number_format(accuracy=0.01))+
  scale_y_discrete(labels=function(x) str_wrap(x,width=35))+
  scale_x_discrete(labels=function(x) gsub("\\s*\\([^\\)]+\\)","",x))+
  xlab("")+
  
  theme(
    panel.grid.minor=element_blank(),
    axis.ticks = element_line(colour="black"),
    axis.line = element_blank(),
    axis.text.x = element_text(size=14, family='sans',angle=30,hjust=1),
    axis.text.y = element_text(size=14, family='sans'),
    axis.title = element_text(size=14, family='sans'),
    text=element_text(size=14, family='sans'),
    legend.text=element_text(size=14, family='sans'),
    plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
    legend.position="right",
    panel.background = element_rect(fill="transparent"),
    panel.border = element_rect(colour="black"),
    legend.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill="transparent", colour=NA)
  )+
  guides(colour=guide_colorbar(expression(Adj.~italic(P)~value),
                               order=1,
                               frame.colour="black",
                               ticks.colour="black", 
                               frame.linewidth=0.5, 
                               ticks.linewidth=0.5, 
                               label.hjust=0, 
                               label.vjust=0.5, 
                               title.hjust=0.5,
                               # title = element_blank(),
                               title.position = "top",
                               title.theme = element_text(family='sans', size=14),
                               label.theme = element_text(family='sans', size=14)))

ggsave(up, filename="up.svg", device=svg, bg="transparent",
       width=7, height=7, units=c("in"))