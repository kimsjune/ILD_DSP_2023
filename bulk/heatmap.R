## Describes code used to generate gene expression heatmaps

## 1: across fibrosis

# show expression of bulk RNA-seq for genes that were identified in fibrosis/foci/central fibrosis regions
efit_AcrossFibrosis <- readRDS("./../nanostr/dge_final/across_fibrosis/efit_across_fibrosis.rds")

# coef 4: fibrosis CHP v FF IPF
# coef 6: central NSIP v FF IPF

# extract top genes
topGenesAcrossFibrosis <-  topTable(efit_AcrossFibrosis, coef=c(4,6), n=50, p.value=0.05, adjust.method="BH", lfc=1)


# IPF M24, M37
# CHP M3, M35, N6
# NSIP M25, M29, M33
# UNC M26, N4, N5







#############
## heatmap of bulk data expression with genes from nanostring
#############

# import y and meta from limma.R
y <- readRDS("./../bulk/dge_final/y.rds")
meta <- readRDS("./../bulk/dge_final/meta.rds")

# subset cpm count table

lcpm_subset_bulk <- cbind(
  edgeR::cpm(y)[,"M24"], #IPF
  edgeR::cpm(y)[,"M34"],
  edgeR::cpm(y)[,"M37"],
  
  
  
  edgeR::cpm(y)[,"M3"], #CHP
  edgeR::cpm(y)[,"M32"],
  edgeR::cpm(y)[,"M35"],
  edgeR::cpm(y)[,"M38"],
  edgeR::cpm(y)[,"M44"],
  
  edgeR::cpm(y)[,"M25"], #NSIP
  edgeR::cpm(y)[,"M29"],
  edgeR::cpm(y)[,"M33"]
  
)
#scale
lcpm_subset_bulk_scale <- cbind(
  id=rownames(lcpm_subset_bulk),
  data.frame(t(scale(t(lcpm_subset_bulk))))
)

# match returns indices in bulk topTags for top 50 genes from Nanostring
lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis <- lcpm_subset_bulk_scale[
  match(rownames(topGenesAcrossFibrosis), rownames(lcpm_subset_bulk_scale)),
]

# row order should be similar to Nanostring heatmap 

# import factor attribute and apply to $id


across_fibrosis_gene_order <- readRDS("./../nanostr/dge_final/across_fibrosis/across_fibrosis_gene_order.rds")
# re-order top genes to match the row dendrogram order
lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis$id <- factor(lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis$id, 
                                                                   # without reversing, the first
                                                                   # element ends up in the bottom row
                                                                   #levels = rev(geneList.order$id))
                                                                   levels = across_fibrosis_gene_order$levels, ordered=T)
# write.table(lcpm_subset_bulk_scale_topGenes$id, file="geneList.dgeTop.txt", quote=F, col.names=F,
#             row.names=F)

# melt() makes the z-score scaled count table into ggplot compatible form
library(reshape2)

lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis_melt <- reshape2::melt(lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis)


# everything in one segment
hm.bulk_with_topGenesAcrossFibrosis <- ggplot(data= lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis_melt, 
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
hm.topGenesBulkAll_in_topGenesAcrossFibrosisAll.legend <- get_legend(hm.bulk_with_topGenesAcrossFibrosis)

# extract heatmap
hm.topGenesBulkAll_in_topGenesAcrossFibrosisAll.noLegend <- hm.bulk_with_topGenesAcrossFibrosis + theme(legend.position="none")


hm.topGenesBulkAll_in_topGenesAcrossFibrosisAll.resized <- set_panel_size(
  hm.topGenesBulkAll_in_topGenesAcrossFibrosisAll.noLegend,
  # width and height of heatmap is determined directly
  # by ncol/nrow in mm, multiplied by a scaling factor of 2
  # reshape() is to unmelt 
  # subtract one column because it's "id" (genes)
  width=2.5*unit(ncol(reshape(lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis_melt, direction="wide", idvar="id", timevar="variable"))-1, "mm"),
  height=2.5*unit(nrow(reshape(lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis_melt, direction="wide", idvar="id", timevar="variable")), "mm"),
  margin=unit(3,"cm"))


ggsave(plot=hm.topGenesBulkAll_in_topGenesAcrossFibrosisAll.resized, 
       filename="hm.topGenesBulkAll_in_topGenesAcrossFibrosisAll.svg", bg="transparent", device=svg,
       # width/height here does not change the ratio;
       # needs to be large enough so that the heatmap doesn't get cut
       width=3, height=6, units=c("in"))


hm.topGenesBulkAll_in_topGenesAcrossFibrosisAll.legend <- as_ggplot(hm.topGenesBulkAll_in_topGenesAcrossFibrosisAll.legend)
ggsave(plot=hm.topGenesBulkAll_in_topGenesAcrossFibrosisAll.legend, 
       # scale is small enough, so width/height doesn't need to be
       filename="hm.topGenesBulkAll_in_topGenesAcrossFibrosisAll.legend.svg", bg="transparent", device=svg,
       width=1, height=2, units=c("in"))

## 2: across hallmark annoations

# show expression of bulk RNA-seq for genes that were identified in hallmark analysis
set.seed(0)
efit_PathDistinctTypes <- readRDS("./../nanostr/dge_final/path_distinct_types/efit_path_distinct_types.rds")

# coef 4: fibrosis CHP v FF IPF
# coef 6: central NSIP v FF IPF

# extract top genes
topGenesPathDistinctTypes <-  topTable(efit_PathDistinctTypes, coef=c(1:3), n=50, p.value=0.05, adjust.method="BH", lfc=1.3)


# IPF M24, M37
# CHP M3, M35, N6
# NSIP M25, M29, M33
# UNC M26, N4, N5







#############
## heatmap of bulk data expression with genes from nanostring
#############

# import y and meta from limma.R
y <- readRDS("./../bulk/dge_final/y.rds")
meta <- readRDS("./../bulk/dge_final/meta.rds")

# subset cpm count table

lcpm_subset_bulk <- cbind(
  edgeR::cpm(y)[,"M24"], #IPF
  edgeR::cpm(y)[,"M34"],
  edgeR::cpm(y)[,"M37"],
  
  
  
  edgeR::cpm(y)[,"M3"], #CHP
  edgeR::cpm(y)[,"M32"],
  edgeR::cpm(y)[,"M35"],
  edgeR::cpm(y)[,"M38"],
  edgeR::cpm(y)[,"M44"],
  
  edgeR::cpm(y)[,"M25"], #NSIP
  edgeR::cpm(y)[,"M29"],
  edgeR::cpm(y)[,"M33"]
)
#scale
lcpm_subset_bulk_scale <- cbind(
  id=rownames(lcpm_subset_bulk),
  data.frame(t(scale(t(lcpm_subset_bulk))))
)

# match returns indices in bulk topTags for top 50 genes from Nanostring
lcpm_subset_bulk_scale_with_topGenesPathDistinctTypes <- lcpm_subset_bulk_scale[
  na.omit(match(rownames(topGenesPathDistinctTypes), rownames(lcpm_subset_bulk_scale))),
]

# row order should be similar to Nanostring heatmap 

# import factor attribute and apply to $id


path_distinct_types_gene_order <- readRDS("./../nanostr/dge_final/path_distinct_types/path_distinct_types_gene_order.rds")
# re-order top genes to match the row dendrogram order
lcpm_subset_bulk_scale_with_topGenesPathDistinctTypes$id <- factor(lcpm_subset_bulk_scale_with_topGenesPathDistinctTypes$id, 
                                                              # without reversing, the first
                                                              # element ends up in the bottom row
                                                              #levels = rev(geneList.order$id))
                                                              levels = path_distinct_types_gene_order$levels, ordered=F)
# write.table(lcpm_subset_bulk_scale_topGenes$id, file="geneList.dgeTop.txt", quote=F, col.names=F,
#             row.names=F)

# melt() makes the z-score scaled count table into ggplot compatible form
library(reshape2)

lcpm_subset_bulk_scale_with_topGenesPathDistinctTypes_melt <- melt(lcpm_subset_bulk_scale_with_topGenesPathDistinctTypes)


# everything in one segment
hm.bulk_with_topGenesPathDistinctTypes <- ggplot(data= lcpm_subset_bulk_scale_with_topGenesPathDistinctTypes_melt, 
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
hm.topGenesBulkAll_in_topGenesPathDistinctTypesAll.legend <- get_legend(hm.bulk_with_topGenesPathDistinctTypes)

# extract heatmap
hm.topGenesBulkAll_in_topGenesPathDistinctTypesAll.noLegend <- hm.bulk_with_topGenesPathDistinctTypes + theme(legend.position="none")


hm.topGenesBulkAll_in_topGenesPathDistinctTypesAll.resized <- set_panel_size(
  hm.topGenesBulkAll_in_topGenesPathDistinctTypesAll.noLegend,
  # width and height of heatmap is determined directly
  # by ncol/nrow in mm, multiplied by a scaling factor of 2
  # reshape() is to unmelt 
  # subtract one column because it's "id" (genes)
  width=2.5*unit(ncol(reshape(lcpm_subset_bulk_scale_with_topGenesPathDistinctTypes_melt, direction="wide", idvar="id", timevar="variable"))-1, "mm"),
  height=2.5*unit(nrow(reshape(lcpm_subset_bulk_scale_with_topGenesPathDistinctTypes_melt, direction="wide", idvar="id", timevar="variable")), "mm"),
  margin=unit(3,"cm"))


ggsave(plot=hm.topGenesBulkAll_in_topGenesPathDistinctTypesAll.resized, 
       filename="hm.topGenesBulkAll_in_topGenesPathDistinctTypesAll.svg", bg="transparent",
       # width/height here does not change the ratio;
       # needs to be large enough so that the heatmap doesn't get cut
       width=3, height=5, units=c("in"))


hm.topGenesBulkAll_in_topGenesPathDistinctTypesAll.legend <- as_ggplot(hm.topGenesBulkAll_in_topGenesPathDistinctTypesAll.legend)
ggsave(plot=hm.topGenesBulkAll_in_topGenesPathDistinctTypesAll.legend, 
       # scale is small enough, so width/height doesn't need to be
       filename="hm.topGenesBulkAll_in_topGenesPathDistinctTypesAll.legend.svg", bg="transparent",
       width=1, height=2, units=c("in"))