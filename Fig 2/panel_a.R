## Depends on standR/limma-voom.R

fit <- readRDS(file= "./../fit.rds")
design <- readRDS(file= "./../design.rds")
dge_all <- readRDS(file="./../dge_all.rds")

con <- makeContrasts(
  fibroblast_v_fibrosis = fibroblast_IPF -fibrosis_IPF,
  fibrosis_v_neutral = fibrosis_IPF - neutral_IPF,
  fibroblast_v_neutral = fibroblast_IPF - neutral_IPF,
  #fibrosis_v_lymphoid = fibrosis_IPF - lymphoid_IPF,
  #lymphoid_v_neutral = lymphoid_IPF - neutral_IPF,
  levels=colnames(design))


### Using batch corrected data with RUV

library(standR)
library(ggplot2)
library(ggpubr)
# import spe_ruv from standR_normalization.R
spe_ruv <- readRDS("./../spe_ruv.rds")

# subset/extract IPF type then do PCA
# https://bioconductor.org/packages/devel/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html

spe_ruv_subset<- spe_ruv[, spe_ruv$type == "IPF"]
spe_ruv_subset<- scater::runPCA(spe_ruv_subset)
pca_ruv_results_subset<- reducedDim(spe_ruv_subset, "PCA")

pca.plot <- drawPCA(spe_ruv_subset, precomputed=pca_ruv_results_subset)+
  #geom_point(colour="black", pch=21, size=2.5)+
  geom_point(aes(shape=anno, fill=anno), size=3, colour="black")+
  scale_shape_manual("Condition: IPF",
    values=c(21:24),
                     breaks=c("neutral","lymphoid","fibrosis","fibroblast"),
                     labels=c("Uninvolved", "Lymphoid", "Fibrosis", "Fibroblastic foci"))+
  scale_fill_manual("Condition: IPF",
        values=c('#FFC2FF',
                 '#FB8EBD',
                 "#C85E78",
                 "#90353B"),
        breaks=c("neutral","lymphoid","fibrosis","fibroblast"),
        labels=c("Uninvolved", "Lymphoid", "Fibrosis", "Fibroblastic foci"))+

      scale_y_continuous(
        labels=scales::number_format(accuracy=0.1))+
      scale_x_continuous(
        labels=scales::number_format(accuracy=0.1))+
      theme_bw()+
      theme(panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            axis.text = element_text(color="black",size=14),
            axis.ticks = element_line(colour="black"),
            axis.title = element_text(size=14),
            #legend.position="none",
            legend.title=element_text(size=14, vjust=0.5, hjust=0.5),
            #legend.title=element_blank(),
            legend.text=element_text(size=14, vjust=0.5, hjust=0),
            plot.margin=unit(c(1,1,1,1),"cm"),
            plot.background=element_rect(fill="transparent", colour=NA),
            panel.border = element_rect(colour="black"),
            panel.background=element_rect(fill="transparent", colour=NA),
            legend.background=element_rect(fill="transparent", colour=NA),
            legend.box.background=element_rect(fill="transparent", colour=NA),
            legend.key=element_rect(fill="transparent", colour=NA),
            legend.position = 'bottom',
            aspect.ratio=1)+
  guides(fill=guide_legend(nrow=2, byrow=T, title.position="top"), shape=guide_legend(nrow=2,byrow=T, title.position="left"))

  

pca.plot_legend <- get_legend(pca.plot)


pca.plot.noLegend <- pca.plot + theme(legend.position="none")

pca.plot.resized <- set_panel_size(
  pca.plot.noLegend,
  width=unit(2.5,"in"), height=unit(2.5,"in")
)
ggsave(plot=pca.plot.resized, 
       filename="PCA_IPF.svg", bg="transparent", device=svg,
       width=5, height=3, units=c("in"))

pca.plot.legend <- as_ggplot(pca.plot_legend)
ggsave(plot=pca.plot.legend, 
       filename="PCA_legend.svg", bg="transparent", device=svg,
       width=5, height=2, units=c("in"))