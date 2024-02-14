## Depends on Fig 3/panel_b.R

# import clustered gene order from across fibrosis comparison
gene_order <- readRDS("./../across_fibrosis/gene_order_xfibrosis.rds")

# extract juts the ones upregulated in central fibrosis compared to rest
gene_order_top <- gene_order[1:(length(gene_order)-12)]

lcpm_subset_annoType_topGenes <- lcpm_subset_annoType[gene_order_top  ,]




# scale, no ID included for heatmap of clusters
lcpm_subset_annoType_topGenes_scale <- t(scale(t(lcpm_subset_annoType_topGenes)))


# scale, include ID for heatmap of genes
lcpm_subset_ID_annoType_topGenes_scale <- cbind(
  id=rownames(lcpm_subset_annoType_topGenes),
  data.frame(t(scale(t(lcpm_subset_annoType_topGenes))))
)


## col dendrogram
col.dendro <- as.dendrogram(hclust(d=dist(t(lcpm_subset_annoType_topGenes_scale)),method="complete"))


# maybe use dendextend instead of ggplot?
dxt.col <- col.dendro %>% 
  #dendextend::set("branches_k_color", value=c("#c85200","#366785"), k=2) %>%
  dendextend::set("labels"," ")  %>%
  dendextend::set("branches_lwd", 0.25) 


dxt.col.gg <- as.ggdend(dxt.col)  

dxt.col.ggPlot <- ggplot(dxt.col.gg,  theme=NULL) +
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    # axis.text.y=element_text(family='sans', face='italic',
    #                          colour='black', size=6, hjust="left"),
    axis.ticks.y=element_blank(),
    axis.title = element_blank(),
    #legend.text=element_text(size=14),
    #legend.title=element_blank(),
    panel.grid = element_blank(),
    plot.margin=unit(c(1,1,1,1),"mm"),
    plot.background=element_rect(fill="transparent", colour=NA),
    panel.background=element_rect(fill="transparent", colour=NA),
    panel.border=element_blank())


ggsave(plot=dxt.col.ggPlot, filename="hm.dge50.col.dendro.svg", bg="transparent", device=svg,
       width=4, height=0.3, units=c("in"))

col.dendro.plot <- ggdendrogram(col.dendro)+
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank())
ggsave(plot=col.dendro.plot,
       filename="hm.dgeAll.col.dendro.svg", bg="transparent", device=svg,
       width=8, height=1, units=c("in"))

order.col <- order.dendrogram(col.dendro)

# define clusters 
col.clusters <- cutree(col.dendro, k=3, order_clusters_as_data=F)


## row dendrogram
# need to do this first to determine the order of genes in lcpm_subset_scale_topGenes
dendro <- as.dendrogram(hclust(d=dist(x=lcpm_subset_annoType_topGenes_scale)))
dendro_plot <- ggdendrogram(data=dendro, rotate=T) +
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank())

dendro.order <- order.dendrogram(dendro)





lcpm_subset_ID_annoType_topGenes_scale$id <- factor(lcpm_subset_ID_annoType_topGenes_scale$id,
                                                    # without reversing, the first
                                                    # element ends up in the bottom row
                                                    #levels = rev(geneList.dge50.order$id))
                                                    levels = lcpm_subset_ID_annoType_topGenes_scale$id[dendro.order], ordered=T)


# melt() makes the z-score scaled count table into ggplot compatible form
library(reshape2)

lcpm_subset_ID_annoType_topGenes_scale_melt <- reshape2::melt(lcpm_subset_ID_annoType_topGenes_scale)

# variables (samples) must be ordered
lcpm_subset_ID_annoType_topGenes_scale_melt$variable <- with(lcpm_subset_ID_annoType_topGenes_scale_melt ,
                                                             factor(variable, 
                                                                    levels=names(col.clusters)))


# everything in one segment
hm.dge50 <- ggplot(data= lcpm_subset_ID_annoType_topGenes_scale_melt, 
                   aes(x=variable,
                       y=id,
                       fill=value))+
  geom_tile(colour="black", linewidth=0.15)+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(family='sans', face='italic',
                                 colour='black', size=6, hjust="left"),
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"),
        plot.background=element_rect(fill="transparent", colour=NA),
        panel.background=element_rect(fill="transparent", colour=NA),
        panel.border=element_blank(),
        legend.background=element_rect(fill="transparent", colour=NA),
        legend.box.background=element_rect(fill="transparent", colour=NA),
        legend.key=element_rect(fill="transparent", colour=NA))+
  
  scale_fill_viridis_c(option="A",
                       guide=guide_colorbar(frame.colour="black",
                                            ticks.colour="black", 
                                            frame.linewidth=0.5, 
                                            ticks.linewidth=0.5, 
                                            label.hjust=0.5, 
                                            label.vjust=0.5, 
                                            label.theme = element_text(
                                              angle=90)))+#,
  scale_y_discrete(position="right", expand=c(0,0))+
  
  #limits=c(-4,4))+
  coord_equal()




# extract legend
hm.dge50.legend <- get_legend(hm.dge50)

# extract heatmap
hm.dge50.noLegend <- hm.dge50 + theme(legend.position="none")


hm.dge50.resized <- set_panel_size(
  hm.dge50.noLegend,
  # width and height of heatmap is determined directly
  # by ncol/nrow in mm, multiplied by a scaling factor of 2
  # reshape() is to unmelt 
  # subtract one column because it's "id" (genes)
  width=2.5*unit(ncol(reshape(lcpm_subset_ID_annoType_topGenes_scale_melt, direction="wide", idvar="id", timevar="variable"))-1, "mm"),
  height=2.5*unit(nrow(reshape(lcpm_subset_ID_annoType_topGenes_scale_melt, direction="wide", idvar="id", timevar="variable")), "mm"),
  margin=unit(3,"cm"))


ggsave(plot=hm.dge50.resized, 
       filename="hm.dge50.svg", bg="transparent", device=svg,
       # width/height here does not change the ratio;
       # needs to be large enough so that the heatmap doesn't get cut
       width=4.2, height=4, units=c("in"))


hm.dge50.legend <- as_ggplot(hm.dge50.legend)
ggsave(plot=hm.dge50.legend, 
       # scale is small enough, so width/height doesn't need to be
       filename="hm.dge50.legend.svg", bg="transparent", device=svg,
       width=1, height=2, units=c("in"))

##########################
### heatmap as column cluster colour bar
##########################


##### Matching columns to patid


# order.col contains the indices of clustered columns of lcpm_subset_annoType_topGenes_scale
# lcpm_subset_annoType does not have original column names, but lcpm_subset does. 
# lcpm_subset is not a SpatialExperiment class because it was made with cbind
# so need to go back to spe_ruv_subset

patid.cluster <- colData(spe_ruv_subset)[colnames(lcpm_subset)[order.col],"patid"]







### separate column colour bars to get separate legends


bar_anno <- rbind(gsub("\\..*","",names(col.clusters)))
colnames(bar_anno) <- names(col.clusters)
rownames(bar_anno) <- c("anno_type")


bar_anno.melt <- melt(bar_anno)

hm.bar_anno <- ggplot(data=bar_anno.melt,
                      aes(x=Var2, y=Var1, fill=value))+
  geom_tile(linewidth=0.15, colour="black")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(family='sans', #face='italic',
                                 colour='black', size=6, hjust="left"),
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"),
        plot.background=element_rect(fill="transparent", colour=NA),
        panel.background=element_rect(fill="transparent", colour=NA),
        panel.border=element_blank(),
        legend.background=element_rect(fill="transparent", colour=NA),
        legend.box.background=element_rect(fill="transparent", colour=NA),
        legend.key=element_rect(fill="transparent", colour=NA))+
  
  scale_fill_manual(values=c(
    "central_NSIP" = "#e37e00",
    
    "fibroblast_IPF" = "#90353B",
    "fibrosis_IPF" = "#C85E78",
    
    "fibrosis_CHP" = "#1a476f",
    
    "fibroblast_UNC" = "#6e8e84",
    "fibrosis_UNC" = "#6e8e84"
    # "1"= "green2",
    # "2" = "blue",
    # "3" = "orange",
    # "M29"=glasbey()[6],
    #   "M33"=glasbey()[7],
    #   "N4"=glasbey()[13],
    #   "M25"=glasbey()[4],
    #   "N5"=glasbey()[14],
    #   "M35"=glasbey()[8],
    #   "N6"=glasbey()[15],
    #   "M37"=glasbey()[9],
    #   "M3"=glasbey()[2],
    #   "M1"=glasbey()[1],
    #   "M24"=glasbey()[3]
    
  ),
  breaks=c("fibrosis_IPF","fibroblast_IPF",
           "central_NSIP",
           "fibrosis_CHP",
           "fibrosis_UNC","fibroblast_UNC"
           
  ),
  labels=c(
    "central_NSIP" = "NSIP",
    "fibroblast_IPF" = "FF IPF",
    "fibrosis_CHP" = "CHP",
    "fibrosis_IPF" = "IPF",
    "fibrosis_UNC" = "UNC.",
    "fibroblast_UNC" = "FF UNC"
    # "1" = "Cluster 1 ",
    # "2" = "Cluster 2 ",
    # "3" = "Cluster 3 "
  ))+
  guides(
    fill=guide_legend(
      #frame.colour="black",
      override.aes=list(linewidth=0.5), # makes the boxes' border thicker
      #ticks.linewidth=0.5,
      label.hjust=0
    ))+
  # guide=guide_colorbar(frame.colour="black",
  #                      ticks.colour="black",
  #                      frame.linewidth=0.5,
  #                      ticks.linewidth=0.5,
  #                      label.hjust=0.5,
  #                      label.vjust=0.5,
  #                      label.theme = element_text(angle=90)))+#,
  scale_y_discrete(position="right", expand=c(0,0),
                   labels=c(#"cluster"="CLUSTER",
                     "anno_type"="Condition"))+
  
  #limits=c(-4,4))+
  coord_equal()

# extract legend
hm.bar_anno.legend <- get_legend(hm.bar_anno)

# extract heatmap
hm.bar_anno.noLegend <- hm.bar_anno + theme(legend.position="none"
)
hm.bar_anno.resized <- set_panel_size(
  hm.bar_anno.noLegend,
  # width and height of heatmap is determined directly
  # by ncol/nrow in mm, multiplied by a scaling factor of 2
  # reshape() is to unmelt
  # subtract one column because it's "id" (genes)
  width=2.5*unit(ncol(reshape(bar_anno.melt, direction="wide", idvar="Var1", timevar="Var2"))-1, "mm"),
  height=2.5*unit(nrow(reshape(bar_anno.melt, direction="wide", idvar="Var1", timevar="Var2")), "mm"),
  margin=unit(3,"cm"))


ggsave(plot=hm.bar_anno.resized,
       filename="hm.bar_anno.svg", bg="transparent", device=svg,
       # width/height here does not change the ratio;
       # needs to be large enough so that the heatmap doesn't get cut
       width=4.5, height=1, units=c("in"))


hm.bar_anno.legend <- as_ggplot(hm.bar_anno.legend)
ggsave(plot=hm.bar_anno.legend,
       # scale is small enough, so width/height doesn't need to be
       filename="hm.bar_anno.legend.svg", bg="transparent", device=svg,
       width=2, height=2, units=c("in"))

#######


bar_patid <- rbind(gsub("\\..*","",patid.cluster))
colnames(bar_patid) <- names(col.clusters)
rownames(bar_patid) <- c("patid")


bar_patid.melt <- melt(bar_patid)

hm.bar_patid <- ggplot(data=bar_patid.melt,
                       aes(x=Var2, y=Var1, fill=value))+
  geom_tile(linewidth=0.15, colour="black")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(family='sans', #face='italic',
                                 colour='black', size=6, hjust="left"),
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"),
        plot.background=element_rect(fill="transparent", colour=NA),
        panel.background=element_rect(fill="transparent", colour=NA),
        panel.border=element_blank(),
        legend.background=element_rect(fill="transparent", colour=NA),
        legend.box.background=element_rect(fill="transparent", colour=NA),
        legend.key=element_rect(fill="transparent", colour=NA))+
  
  scale_fill_manual(values=c(
    # "central_NSIP" = "#9F7468",
    # 
    # "fibroblast_IPF" = "#90353B",
    # "fibrosis_IPF" = "#C85E78",
    # 
    # "fibrosis_CHP" = "#334B49",
    # 
    # "fibroblast_UNC" = "grey50",
    # "fibrosis_UNC" = "grey10"
    # "1"= "green2",
    # "2" = "blue",
    # "3" = "orange",
    # "M29"=#glasbey()[6],
    #  "M33"=#glasbey()[7],
    "N4"=glasbey()[13],
    # "M25"=#glasbey()[4],
    "M26"=glasbey()[5],
    "N5"=glasbey()[14]
    #  "M35"=#glasbey()[8],
    #  "N6"=#glasbey()[15],
    #  "M37"=#glasbey()[9],
    #  "M3"=#glasbey()[2],
    #  "M1"=#glasbey()[1],
    #  "M24"=#glasbey()[3]
    
    #),
    # labels=c(
    #   "central_NSIP" = "NSIP",
    #   "fibroblast_IPF" = "IPF",
    #   "fibrosis_CHP" = "CHP",
    #   "fibrosis_IPF" = "IPF",
    #   "fibrosis_UNC" = "fibrosis Unc.",
    #   "fibroblast_UNC" = "Unc."
    #   # "1" = "Cluster 1 ",
    #   # "2" = "Cluster 2 ",
    #   # "3" = "Cluster 3 "
  ),
  na.value = "white")+
  guides(
    fill=guide_legend(
      #frame.colour="black",
      override.aes=list(linewidth=0.5), # makes the boxes' border thicker
      #ticks.linewidth=0.5,
      label.hjust=0,
      ncol=1
    ))+
  # guide=guide_colorbar(frame.colour="black",
  #                      ticks.colour="black",
  #                      frame.linewidth=0.5,
  #                      ticks.linewidth=0.5,
  #                      label.hjust=0.5,
  #                      label.vjust=0.5,
  #                      label.theme = element_text(angle=90)))+#,
  scale_y_discrete(position="right", expand=c(0,0),
                   labels=c(#"cluster"="CLUSTER",
                     "patid"="Patient ID"))+
  
  #limits=c(-4,4))+
  coord_equal()

# extract legend
hm.bar_patid.legend <- get_legend(hm.bar_patid)

# extract heatmap
hm.bar_patid.noLegend <- hm.bar_patid + theme(legend.position="none"
)
hm.bar_patid.resized <- set_panel_size(
  hm.bar_patid.noLegend,
  # width and height of heatmap is determined directly
  # by ncol/nrow in mm, multiplied by a scaling factor of 2
  # reshape() is to unmelt
  # subtract one column because it's "id" (genes)
  width=2.5*unit(ncol(reshape(bar_patid.melt, direction="wide", idvar="Var1", timevar="Var2"))-1, "mm"),
  height=2.5*unit(nrow(reshape(bar_patid.melt, direction="wide", idvar="Var1", timevar="Var2")), "mm"),
  margin=unit(3,"cm"))


ggsave(plot=hm.bar_patid.resized,
       filename="hm.bar_patid.svg", bg="transparent", device=svg,
       # width/height here does not change the ratio;
       # needs to be large enough so that the heatmap doesn't get cut
       width=4.5, height=1, units=c("in"))


hm.bar_patid.legend <- as_ggplot(hm.bar_patid.legend)
ggsave(plot=hm.bar_patid.legend,
       # scale is small enough, so width/height doesn't need to be
       filename="hm.bar_patid.legend.svg", bg="transparent", device=svg,
       width=2, height=2, units=c("in"))