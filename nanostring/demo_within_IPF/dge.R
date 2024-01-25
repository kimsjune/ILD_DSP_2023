### following limma-voom pipeline to find DGEs between annotations
set.seed(0)
fit_contrast <- contrasts.fit(fit, contrasts = con)
efit <- eBayes(fit_contrast, robust = TRUE)

# save efit to compare against Blumhagen data
saveRDS(efit, "efit_within_IPF.rds")

results_efit<- decideTests(efit, p.value = 0.05, lfc=1)
summary_efit <- summary(results_efit)

resultsLogFC_efit<- decideTests(efit, p.value = 0.05, lfc=1.5)
summaryLogFC_efit <- summary(resultsLogFC_efit)

topGenes <- topTable(efit, coef=c(1:3), n=50, p.value=0.05, adjust.method="BH", lfc=1, sort.by="F")

topGenesTotal <- topTable(efit, coef=c(1:3), n=200, p.value=0.05, adjust.method="BH")


## RUV4-normalized log CPM whose columns are grouped by annotation
## also 

lcpm_subset <- cbind(
 # id=rownames(assay(spe_ruv_subset,2)),
    assay(spe_ruv_subset,2)[, colData(spe_ruv_subset)$anno=="neutral"], # uninvolved regions
    assay(spe_ruv_subset,2)[, colData(spe_ruv_subset)$anno=="fibrosis"],
    assay(spe_ruv_subset,2)[, colData(spe_ruv_subset)$anno=="fibroblast"])

# scale
lcpm_subset_scale <- cbind(
  id=rownames(lcpm_subset),
  data.frame(t(scale(t(lcpm_subset))))
)

# subset topGenes only

lcpm_subset_scale_topGenes <- lcpm_subset_scale[rownames(topGenes),]

lcpm_subset_scale_topGenesTotal <- lcpm_subset_scale[rownames(topGenesTotal),]



## row dendrogram
# need to do this first to determine the order of genes in lcpm_subset_scale_topGenes

dendro <- as.dendrogram(hclust(d=dist(x=lcpm_subset_scale_topGenes)))
# dendro_plot <- ggdendrogram(data=dendro, rotate=T) +
#   theme(axis.text.y=element_blank(),
#         axis.text.x=element_blank())

# maybe use dendextend instead of ggplot?
dxt <- dendro %>% 
  dendextend::set("labels","")  %>%
  dendextend::set("branches_lwd", 0.25) %>%
  dendextend::set("branches_k_color",k=2, value=c("#c85200","#366785"))# %>%
  #rotate(50:1)
  
dxt.gg <- as.ggdend(dxt)  

dxt.ggPlot <- ggplot(dxt.gg,  theme=NULL) +
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
         panel.border=element_blank())+
   coord_flip()

# ggsave(plot=dendro_plot, 
#        # scale is small enough, so width/height doesn't need to be
#        filename="hm.dge50.dendro.svg", bg="transparent",
#        width=1, height=6, units=c("in"))

ggsave(plot=dxt.ggPlot, filename="hm.dge50.dendroGG.svg", bg="transparent", device=svg,
       width=0.3, height=4, units=c("in"))

dendro.order <- order.dendrogram(dendro)

# re-order top genes to match the row dendrogram order
lcpm_subset_scale_topGenes$id <- factor(lcpm_subset_scale_topGenes$id, 
                                                  # without reversing, the first
                                                  # element ends up in the bottom row
                                                  #levels = rev(geneList.dge50.order$id))
                                   levels = lcpm_subset_scale_topGenes$id[dendro.order], ordered=T)
# write.table(lcpm_subset_scale_topGenes$id, file="geneList.dgeTop.txt", quote=F, col.names=F,
#             row.names=F)

# melt() makes the z-score scaled count table into ggplot compatible form
library(reshape2)

lcpm_subset_scale_topGenes_melt <- melt(lcpm_subset_scale_topGenes)


# everything in one segment
hm.dge50 <- ggplot(data= lcpm_subset_scale_topGenes_melt, 
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
  width=2.5*unit(ncol(reshape(lcpm_subset_scale_topGenes_melt, direction="wide", idvar="id", timevar="variable"))-1, "mm"),
  height=2.5*unit(nrow(reshape(lcpm_subset_scale_topGenes_melt, direction="wide", idvar="id", timevar="variable")), "mm"),
  margin=unit(3,"cm"))


ggsave(plot=hm.dge50.resized, 
       filename="hm.dge50.svg", bg="transparent", device=svg,
       # width/height here does not change the ratio;
       # needs to be large enough so that the heatmap doesn't get cut
       width=4, height=6, units=c("in"))


hm.dge50.legend <- as_ggplot(hm.dge50.legend)
ggsave(plot=hm.dge50.legend, 
       # scale is small enough, so width/height doesn't need to be
       filename="hm.dge50.legend.svg", bg="transparent",
       width=1, height=2, units=c("in"))


###################

###
### Single coefficient topTables (for logFC heatmap, Volcano plots, ORA)
###

topGenes_fibroblast_v_fibrosis <- topTable(efit, coef=1, n=Inf, p.value=0.05, adjust.method="BH", lfc=1)
topGenes_fibrosis_v_neutral <- topTable(efit, coef=2, n=Inf, p.value=0.05, adjust.method="BH", lfc=1)
topGenes_fibroblast_v_neutral <- topTable(efit, coef=3, n=Inf, p.value=0.05, adjust.method="BH", lfc=1)


## In the clustered rows/genes, where do significantly differential genes occur for each coefficient? Do they belong in a particular cluster?

# match single coefficient DGEs to clustered colnames/genes and return boolean for each coef then concatenate by rows 
# coef_df <-data.frame(
#   genes = rep(lcpm_subset_scale_topGenes$id[dendro.order],3),
#   coef=factor(c(rep("coef_fibroblast_v_fibrosis",50),
#          rep("coef_fibrosis_v_neutral",50),
#          rep("coef_fibroblast_v_neutral",50)),
#          levels=c("coef_fibroblast_v_fibrosis",
#                   "coef_fibrosis_v_neutral",
#                   "coef_fibroblast_v_neutral")),
#   bool =c(as.character(lcpm_subset_scale_topGenes$id[dendro.order] %in% rownames(topGenes_fibroblast_v_fibrosis)),
#           as.character(lcpm_subset_scale_topGenes$id[dendro.order] %in% rownames(topGenes_fibrosis_v_neutral)),
#           as.character(lcpm_subset_scale_topGenes$id[dendro.order] %in% rownames(topGenes_fibroblast_v_neutral))))
          

coef_df <- data.frame(
  genes = lcpm_subset_scale_topGenes$id[dendro.order],
  coef_fibroblast_v_fibrosis = topGenes_fibroblast_v_fibrosis[match(lcpm_subset_scale_topGenes$id[dendro.order],rownames(topGenes_fibroblast_v_fibrosis)),"logFC"],
  coef_fibrosis_v_neutral = topGenes_fibrosis_v_neutral[match(lcpm_subset_scale_topGenes$id[dendro.order],rownames(topGenes_fibrosis_v_neutral)),"logFC"],
  coef_fibroblast_v_neutral = topGenes_fibroblast_v_neutral[match(lcpm_subset_scale_topGenes$id[dendro.order],rownames(topGenes_fibroblast_v_neutral)),"logFC"]
)

coef_df_order <- coef_df[lcpm_subset_scale_topGenes$id[dendro.order],]

coef_df_order_melt <- melt(coef_df_order)





####



hm.logFC <- ggplot(data= coef_df_order_melt, 
       aes(x=variable,
           y=genes,
           fill=value
          ))+
  geom_tile(colour="black", linewidth=0.15)+
  # geom_tile(data= ~ subset(.x, is.na(value)),
  #           values=c("NA"="black"))+
  # geom_text(data= ~ subset(.x, is.na(value)), aes(label="|"),
  #           angle=-45, 
  #           family="sans Nova Light",
  #           size=3,
  #           nudge_x=0.1,
  #           nudge_y=0.1
  #           )+
  # geom_segment(data= ~ subset(.x, is.na(value)),
  #              aes(x=xmin, xend=xmax, y=ymin, yend=ymax), 
  #              color="black", linewidth=2)+
  geom_tile_pattern(data = ~ subset(.x, is.na(value)),
                    show.legend=TRUE,
                    pattern_color=NA,
                    pattern_fill="black",
                    pattern_angle=45,
                    pattern="stripe",
                    pattern_size=1.2,
                    pattern_spacing=0.05)+
  geom_tile(colour="black", linewidth=0.15)+
  
  theme_bw()+
  theme(axis.text.x.top=element_text(size=6, colour="black",vjust=1, hjust=0.5, margin=margin(0,0,-3,0) ),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        # axis.text.y=element_text(family='sans', face='italic',
        #                          colour='black', size=6, hjust="left"),
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        legend.text=element_text(size=14),
        legend.title=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"),
        plot.background=element_rect(fill="transparent", colour=NA),
        panel.background=element_rect(fill="transparent", colour=NA),
        panel.grid = element_blank(),
        panel.border=element_blank(),
        legend.background=element_rect(fill="transparent", colour=NA),
        legend.box.background=element_rect(fill="transparent", colour=NA),
        legend.box="vertical",
        legend.key=element_rect(fill="transparent", colour=NA))+
  # theme(#axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank(),
  #       axis.text.y=element_text(family='sans', face='italic',
  #                                colour='black', size=6, hjust="left"),
  #       axis.ticks.y=element_blank(),
  #       axis.title = element_blank(),
  #       legend.text=element_text(size=14),
  #       legend.title=element_blank(),
  #       plot.margin=unit(c(1,1,1,1),"cm"),
  #       plot.background=element_rect(fill="transparent", colour=NA),
  #       panel.background=element_rect(fill="transparent", colour=NA),
  #       panel.border=element_rect(colour="black"),
  #       legend.background=element_rect(fill="transparent", colour=NA),
  #       legend.box.background=element_rect(fill="transparent", colour=NA),
  #       legend.key=element_rect(fill="transparent", colour=NA))+
  scale_fill_scico(
                   palette="vikO", na.value="transparent", 
                   limits=c(-6,6),
                   guide=guide_colorbar(
                                        frame.colour="black",
                                        ticks.colour="black", 
                                        frame.linewidth=0.5, 
                                        ticks.linewidth=0.5, 
                                        label.hjust=0.5, 
                                        label.vjust=0.5, 
                                        label.theme = element_text(
                                          angle=90)))+
  scale_x_discrete(labels=c("A","B","C"),
                   position="top")
    # guide=guide_legend(frame.colour="black",
    #                    # ticks.colour="black", 
    #                     frame.linewidth=0.5, 
    #                    # ticks.linewidth=0.5, 
    #                     label.hjust=0.5, 
    #                     label.vjust=0.5, 
    #                     label.theme = element_text(
    #                       angle=90))
  
  
  # scale_fill_viridis_c(option="A",
  #                      guide=guide_colorbar(frame.colour="black",
  #                                           ticks.colour="black", 
  #                                           frame.linewidth=0.5, 
  #                                           ticks.linewidth=0.5, 
  #                                           label.hjust=0.5, 
  #                                           label.vjust=0.5, 
  #                                           label.theme = element_text(
  #                                             angle=90)))+#,
  #scale_y_discrete(position="right", expand=c(0,0))+
  
  #limits=c(-4,4))+
  #coord_equal()


# extract legend
hm.logFC.legend <- get_legend(hm.logFC)

# extract heatmap
hm.logFC.noLegend <- hm.logFC + theme(legend.position="none")


hm.logFC.resized <- set_panel_size(
  hm.logFC.noLegend,
  # width and height of heatmap is determined directly
  # by ncol/nrow in mm, multiplied by a scaling factor of 2
  # reshape() is to unmelt 
  # subtract one column because it's "id" (genes)
  width=2.5*unit(ncol(reshape(coef_df_order_melt, direction="wide", idvar="genes", timevar="variable"))-1, "mm"),
  height=2.5*unit(nrow(reshape(coef_df_order_melt, direction="wide", idvar="genes", timevar="variable")), "mm"),
  margin=unit(3,"cm"))


ggsave(plot=hm.logFC.resized, 
       filename="hm.logFC.svg", bg="transparent", device=svg,
       # width/height here does not change the ratio;
       # needs to be large enough so that the heatmap doesn't get cut
       width=1, height=6, units=c("in"))


hm.logFC.legend <- as_ggplot(hm.logFC.legend)
ggsave(plot=hm.logFC.legend, 
       # scale is small enough, so width/height doesn't need to be
       filename="hm.logFC.legend.svg", bg="transparent", device=svg,
       width=1, height=2, units=c("in"))

############
### patid as colour bar
############

# colnames(lcpm_subset) returns full slide names, which are rownames of colData(spe_ruv_subset)
# then select patid column of colData

patid_subset <- cbind(
  val=colData(spe_ruv_subset)[colnames(lcpm_subset),"patid"])

# x axis must NOT be numerical because ggplot treats it as continuous and width of tiles won't line up with heatmaps...
rownames(patid_subset) <- paste0(colData(spe_ruv_subset)[colnames(lcpm_subset),"patid"],".",c(1:length(colData(spe_ruv_subset)[colnames(lcpm_subset),"patid"])))
  

patid_subset_melt <- reshape2::melt(t(patid_subset))

colbar <- ggplot(data=patid_subset_melt,
                 aes(y=Var1, x=Var2, fill=value))+
  geom_tile(aes(fill=value), linewidth=0.15, colour="black")+
  theme_bw()+
  # theme(axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank(),
  #       axis.text.y=element_text(family='sans', #face='italic',
  #                                colour='black', size=6, hjust="left"),
  #       
  #       axis.ticks.y=element_blank(),
  #       axis.title = element_blank(),
  #       legend.text=element_text(size=14),
  #       legend.title=element_blank(),
  #       plot.margin=unit(c(1,1,1,1),"cm"),
  #       plot.background=element_rect(fill="transparent", colour=NA),
  #       panel.background=element_rect(fill="transparent", colour=NA),
  #       panel.border=element_blank(),
  #       panel.grid = element_blank(),
  #       legend.background=element_rect(fill="transparent", colour=NA),
  #       legend.box.background=element_rect(fill="transparent", colour=NA),
  #       legend.key=element_rect(fill="transparent", colour=NA))+
theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_text(family='sans', face='plain',
                               colour='black', size=6),
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
    "M1" = "#0000FF",
    "M37" = "#9A4D42",
    "M24" = "#00FF00"

  ),
  labels=c("val"="Patient ID"

  ))+
  guides(
    fill=guide_legend(
      #frame.colour="black",
      override.aes=list(linewidth=0.5), # makes the boxes' border thicker
      #ticks.linewidth=0.5,
      label.hjust=0,
      nrow=1, byrow=T
    ))+
  # guide=guide_colorbar(frame.colour="black",
  #                      ticks.colour="black", 
  #                      frame.linewidth=0.5, 
  #                      ticks.linewidth=0.5, 
  #                      label.hjust=0.5, 
  #                      label.vjust=0.5, 
  #                      label.theme = element_text(angle=90)))+#,
  scale_y_discrete( expand=c(0,0), # appreas on the left
                   labels=c("val"="Patient ID"),
                   position="right")

# extract legend
colbar.legend <- get_legend(colbar)

# extract heatmap
colbar.noLegend <- colbar + theme(legend.position="none"
)
colbar.resized <- set_panel_size(
  colbar.noLegend,
  # width and height of heatmap is determined directly
  # by ncol/nrow in mm, multiplied by a scaling factor of 2
  # reshape() is to unmelt 
  # subtract one column because it's "id" (genes)
  width=2.5*unit(ncol(reshape(patid_subset_melt, direction="wide", idvar="Var1", timevar="Var2"))-1, "mm"), # no -1 because there is no id?
  height=2.5*unit(nrow(reshape(patid_subset_melt, direction="wide", idvar="Var1", timevar="Var2")), "mm"),
  margin=unit(3,"cm"))


ggsave(plot=colbar.resized, 
       filename="colbar.svg", bg="transparent", device=svg,
       # width/height here does not change the ratio;
       # needs to be large enough so that the heatmap doesn't get cut
       width=3.5, height=0.5, units=c("in"))
ggsave(plot=colbar.legend, 
       filename="colbarLegend.svg", bg="transparent", device=svg,
       # width/height here does not change the ratio;
       # needs to be large enough so that the heatmap doesn't get cut
       width=3, height=1, units=c("in"))

