### custom contrasts
# need to average two or more annotations together into one
con_custom <- makeContrasts(
  fibroblast_IPF_v_rest = fibroblast_IPF - (neutral_IPF+
                                              
                                              fibrosis_IPF)/2,
  neutral_IPF_v_rest = neutral_IPF - (fibrosis_IPF+
                                          
                                          fibroblast_IPF)/2,
  
  
  levels=design)


# make contrasts
set.seed(0)
fit_contrast_custom <- contrasts.fit(fit, contrasts = con_custom)
efit_custom <- eBayes(fit_contrast_custom, robust = TRUE)


# save genes and log FC
topGenes_turquoise <- topTable(efit_custom, coef=1, n=Inf)
topGenes_orange <- topTable(efit_custom, coef=2, n=Inf)

# need a named list of entrez gene id and logFC
geneList_turquoise <- as.numeric(topGenes_turquoise[,'logFC'])
names(geneList_turquoise) <- unname(mapIds(org.Hs.eg.db, rownames(topGenes_turquoise), 'ENTREZID', 'SYMBOL'))
geneList_orange <- as.numeric(topGenes_orange[,'logFC'])
names(geneList_orange) <- unname(mapIds(org.Hs.eg.db, rownames(topGenes_orange), 'ENTREZID', 'SYMBOL'))



# get the gene names of each cluster then convert to entrez
orange <- names(cutree(dendro, k=2, order_clusters_as_data=F)[cutree(dendro, k=2, order_clusters_as_data=F)==1])
turquoise <-names(cutree(dendro, k=2, order_clusters_as_data=F)[cutree(dendro, k=2, order_clusters_as_data=F)==2])  

orange_entrez <- mapIds(org.Hs.eg.db, orange, 'ENTREZID', 'SYMBOL')
turquoise_entrez <- mapIds(org.Hs.eg.db, turquoise, 'ENTREZID', 'SYMBOL')


library("org.Hs.eg.db")

library(enrichplot)
library(DOSE)
library(clusterProfiler)
ego_orange <- enrichGO(orange_entrez, 
                       OrgDb=org.Hs.eg.db,
                       ont="BP")
egox_orange <- setReadable(ego_orange, 'org.Hs.eg.db', 'ENTREZID')

ego_turquoise <- enrichGO(turquoise_entrez, 
                          OrgDb=org.Hs.eg.db,
                          ont="BP")
egox_turquoise <- setReadable(ego_turquoise, 'org.Hs.eg.db', 'ENTREZID')

############ cnetplot
# turquoise
cnetplot_turquoise <-cnetplot(egox_turquoise, foldChange=sort(geneList_turquoise, decreasing=T), 
                              colorEdge = F,
                              circular=F, layout='graphopt')+
  scale_colour_scico(palette="vikO", na.value="transparent", limits=c(-3,3))+
  
  theme(
    
    text=element_text(size=20, family='sans'),
    legend.text=element_text(size=14, family='sans'),
    plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")
  )+
  guides(colour=guide_colorbar(expression(log[2]~fold~change),
                               frame.colour="black",
                               ticks.colour="black", 
                               frame.linewidth=0.5, 
                               ticks.linewidth=0.5, 
                               label.hjust=0.5, 
                               label.vjust=0.5, 
                               title.hjust=0.5,
                               # title = element_blank(),
                               title.position = "right",
                               title.theme = element_text(family='sans', size=14, angle=90),
                               label.theme = element_text(family='sans', size=14, angle=90)),
         size=guide_legend("Size",
                           #title=element_blank(),
                           label.hjust=0.5,
                           label.vjust=0.5,
                           title.hjust=0.5,
                           title.position = "right",
                           title.theme = element_text(family='sans', size=14, angle=90),
                           label.theme = element_text(family='sans', size=14, angle=90)),
         legend=guide_bins(title=element_blank()))

cnetplot_turquoise_legend <- get_legend(cnetplot_turquoise)

cnetplot_turquoise.noLegend <- cnetplot_turquoise+ theme(legend.position="none")

cnetplot_turquoise.resized <- set_panel_size(
  cnetplot_turquoise.noLegend,
  width=unit(5,"in"), height=unit(6,"in")
)
ggsave(plot=cnetplot_turquoise.resized, 
       filename="cnetplot_turquoise.resized_graphopt.svg", bg="transparent", device=svg,
       width=5, height=6, units=c("in"))

ggsave(plot=cnetplot_turquoise_legend,
       filename="cnetplot_turquoise_legend.svg", bg="transparent",device=svg,
       width=1.2, height=4, units=c("in"))


# orange
cnetplot_orange <-cnetplot(egox_orange, foldChange=sort(geneList_orange, decreasing=T), 
                           colorEdge = F,
                           circular=F, layout='graphopt')+
  scale_colour_scico(palette="vikO", na.value="transparent", limits=c(-3,3))+
  
  theme(
    
    text=element_text(size=20, family='sans'),
    legend.text=element_text(size=14, family='sans'),
    plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")
  )+
  guides(colour=guide_colorbar(expression(log[2]~fold~change),
                               order=2,
                               frame.colour="black",
                               ticks.colour="black", 
                               frame.linewidth=0.5, 
                               ticks.linewidth=0.5, 
                               label.hjust=0.5, 
                               label.vjust=0.5, 
                               title.hjust=0.5,
                               # title = element_blank(),
                               title.position = "right",
                               title.theme = element_text(family='sans', size=14, angle=90),
                               label.theme = element_text(family='sans', size=14, angle=90)),
         size=guide_legend("Size",
                           #title=element_blank(),
                           order=1,
                           label.hjust=0.5,
                           label.vjust=0.5,
                           title.hjust=0.5,
                           title.position = "right",
                           title.theme = element_text(family='sans', size=14, angle=90),
                           label.theme = element_text(family='sans', size=14, angle=90)),
         legend=guide_bins(title=element_blank()))

cnetplot_orange_legend <- get_legend(cnetplot_orange)

cnetplot_orange.noLegend <- cnetplot_orange+ theme(legend.position="none")

cnetplot_orange.resized <- set_panel_size(
  cnetplot_orange.noLegend,
  width=unit(5,"in"), height=unit(4,"in")
)
ggsave(plot=cnetplot_orange.resized, 
       filename="cnetplot_orange.resized_graphopt.svg", bg="transparent", device=svg,
       width=5, height=4, units=c("in"))

ggsave(plot=cnetplot_orange_legend,
       filename="cnetplot_orange_legend.svg", bg="transparent",device=svg,
       width=1.2, height=4, units=c("in"))









##########
# dotplot 
# turquoise
dot_turquoise <- dotplot(ego_turquoise, showCategory=5)+
  scale_colour_viridis(limits=c(as.double(1e-20), 5e-2),
                       trans="log",
                       breaks=c(1e-20, 1e-10, 1e-5, 5e-2))+
  scale_size_continuous(limits=c(1,15),
                        breaks=c(5,10,15))+
  
  scale_x_continuous(
    labels=scales::number_format(accuracy=0.01))+
  
  theme(
    panel.grid.minor=element_blank(),
    axis.ticks = element_line(colour="black"),
    axis.line = element_blank(),
    axis.text.x = element_text(size=14, family='sans', angle=30,hjust=1),
    axis.text.y = element_text(size=14, family='sans'),
    axis.title = element_text(size=14, family='sans'),
    text=element_text(size=14, family='sans'),
    legend.text=element_text(size=14, family='sans'),
    plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
    legend.position="bottom",
    panel.background = element_rect(fill="transparent"),
    panel.border = element_rect(colour="black"),
    legend.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill="transparent", colour=NA)
  )+
  guides(colour=guide_colorbar("Adj. P val",
                               order=1,
                               frame.colour="black",
                               ticks.colour="black", 
                               frame.linewidth=0.5, 
                               ticks.linewidth=0.5, 
                               label.hjust=1, 
                               label.vjust=1, 
                               title.hjust=0.5,
                               # title = element_blank(),
                               title.position = "top",
                               title.theme = element_text(family='sans', size=14),
                               label.theme = element_text(family='sans', size=14, angle=45)),
         size=guide_legend("Count",
                           #title=element_blank(),
                           label.hjust=0.5,
                           label.vjust=0.5,
                           title.hjust=0.5,
                           title.position = "top",
                           title.theme = element_text(family='sans', size=14), 
                           label.theme = element_text(family='sans', size=14)))

dot_turquoise.noLegend <- dot_turquoise  + theme(legend.position="none")
dot_turquoise.resized <- set_panel_size(
  dot_turquoise.noLegend,
  width=unit(2,"in"), height=unit(2.5,"in")
)
ggsave(plot=dot_turquoise.resized, 
       filename="dot_turquoise.resized.svg", bg="transparent", device=svg,
       width=5, height=4, units=c("in"))

dot_turquoise.legend <- get_legend(dot_turquoise)
dot_turquoise.legend <- as_ggplot(dot_turquoise.legend)
ggsave(plot=dot_turquoise.legend,
       filename="dot_turquoise.legend.svg", bg="transparent",
       device=svg,
       width=6, height=1.5, units=c("in"))


# orange

dot_orange <- dotplot(ego_orange, showCategory=5)+
  scale_colour_viridis(limits=c(as.double(1e-20), 5e-2),
                       trans="log",
                       breaks=c(1e-20, 1e-10, 1e-5, 5e-2))+
  scale_size_continuous(limits=c(1,15),
                        breaks=c(5,10,15))+
  
  scale_x_continuous(
    labels=scales::number_format(accuracy=0.01))+
  
  theme(
    panel.grid.minor=element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour="black"),
    axis.text.x = element_text(size=14, family='sans',angle=30, hjust=1),
    axis.text.y = element_text(size=14, family='sans'),
    axis.title = element_text(size=14, family='sans'),
    text=element_text(size=14, family='sans'),
    legend.text=element_text(size=14, family='sans'),
    plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
    legend.position="bottom",
    panel.background = element_rect(fill="transparent"),
    panel.border = element_rect(colour="black"),
    legend.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill="transparent", colour = "transparent")
  )+
  guides(colour=guide_colorbar("Adj. P val",
                               order=1,
                               frame.colour="black",
                               ticks.colour="black", 
                               frame.linewidth=0.5, 
                               ticks.linewidth=0.5, 
                               label.hjust=1, 
                               label.vjust=1, 
                               title.hjust=0.5,
                               # title = element_blank(),
                               title.position = "top",
                               title.theme = element_text(family='sans', size=14),
                               label.theme = element_text(family='sans', size=14, angle=45)),
         size=guide_legend("Count",
                           #title=element_blank(),
                           label.hjust=0.5,
                           label.vjust=0.5,
                           title.hjust=0.5,
                           title.position = "top",
                           title.theme = element_text(family='sans', size=14), 
                           label.theme = element_text(family='sans', size=14)))

dot_orange.noLegend <- dot_orange  + theme(legend.position="none")
dot_orange.resized <- set_panel_size(
  dot_orange.noLegend,
  width=unit(2,"in"), height=unit(2.5,"in")
)
ggsave(plot=dot_orange.resized, 
       filename="dot_orange.resized.svg", bg="transparent", device=svg,
       width=5, height=4, units=c("in"))
dot_orange.legend <- get_legend(dot_orange)
dot_orange.legend <- as_ggplot(dot_orange.legend)
ggsave(plot=dot_orange.legend,
       filename="dot_orange.legend.svg", bg="transparent",
       device=svg,
       width=6, height=1.5, units=c("in"))

