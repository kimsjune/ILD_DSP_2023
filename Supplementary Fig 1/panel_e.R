### batch correction with RUV

speRUV <- findNCGs(speFinal, batch_name="patid", top_n=200)
for(i in seq(4)){
  spe_ruv <- geomxBatchCorrection(speRUV, factors = "type", 
                                  NCGs = metadata(speRUV)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = patid, title = paste0("k = ", i)))
  
}


spe_ruv <- geomxBatchCorrection(speRUV, factors = "anno_type", 
                                NCGs = metadata(speRUV)$NCGs, k = 4)


saveRDS(spe_ruv, "spe_ruv.rds")


## Relative log2 expression


# make type into a factor
colData(spe_ruv)$type <- factor(colData(spe_ruv)$type, 
                               levels=c("IPF","NSIP","CHP","UNC","NOR"),
                               ordered=T)
# returns colData where rows are sorted by type from IPF to NOR
colData(spe_ruv)[do.call(order, colData(spe_ruv)["type"]),]

RLExpr_sizeNormRuv <- plotRLExpr(
  spe_ruv[,rownames(colData(spe_ruv)[do.call(order, colData(spe_ruv)["type"]),])], 
                                 assay=2, color=type, ordannots="type") +
  geom_boxplot()+
  geom_point(colour="black",size=1)+
  
  scale_colour_manual(
    values=c(
      "CHP"="#1a476f",
      "IPF"= "#90353b",
      "NOR"="#55752f",
      "NSIP"="#e37e00",
      "UNC" ="#6e8e84"
    ),
    breaks=c("IPF","NSIP","CHP","UNC","NOR"))+
  xlab("ROIs")+
  theme(#panel.background = element_rect(fill="grey90"),
    axis.ticks = element_line(),
    text=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    legend.title=element_blank() ,#element_text(size=14, colour="black", face='plain'),
    plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm"))+
  theme(aspect.ratio=1)

ggsave(plot= RLExpr_sizeNormRuv, file="RLExpr_sizeNormRuv.svg", bg="transparent", device=svg,
       width=5, height=4, units=c("in"))



## PCA
spe_ruv <- scater::runPCA(spe_ruv)
pca_results_ruv <- reducedDim(spe_ruv, "PCA")

pca_plot_ruv <- drawPCA(spe_ruv, precomputed=pca_results_ruv)+
  # scale_color_manual(values=c(stata_pal("s2color")(16)))
  # scale_color_tableau(
  #   palette="Classic 20",
  #   type=c("ordered-diverging"),
  #   direction=1
  geom_point(aes(fill=factor(patid)), colour="black", pch=21, size=2.5)+
  scale_fill_manual(guide="legend",
                      name="Patient ID",
                      values=c(unname(glasbey())),
                      breaks=c(mixedsort(unique(spe_q3$patid))))+
  # scale__discrete(breaks=c(mixedsort(unique(spe_q3$patid))))+
  theme(#panel.background = element_rect(fill="grey90"),
    axis.ticks = element_line(),
    text=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    legend.title=element_blank() ,#element_text(size=14, colour="black", face='plain'),
    plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm"))+
  theme(aspect.ratio=1)
ggsave(plot = pca_plot_ruv, file="pca_plot_ruv.svg", device=svg,
       bg="transparent", width=5, height=4, units=c("in"))


## by type 
pca_plot_ruv_type <- drawPCA(spe_ruv, precomputed = pca_results_ruv) +
  # scale_color_manual(values=c(stata_pal("s2color")(16)))+
  scale_fill_manual(
    values=c(
      "CHP"="#1a476f",
      "IPF"= "#90353b",
      "NOR"="#55752f",
      "NSIP"="#e37e00",
      "UNC" ="#6e8e84"
    ),
    breaks=c("IPF","NSIP","CHP","UNC","NOR"))+
  
  geom_point(aes(fill=factor(type)),colour="black", pch=21, size=2.5)+
  # scale_colour_manual(guide="legend",
  #                     name="Patient ID",
  #                     values=c(unname(glasbey())),
  #                     breaks=c(mixedsort(unique(spe_q3$patid))))+
  # scale__discrete(breaks=c(mixedsort(unique(spe_q3$patid))))+
  theme(#panel.background = element_rect(fill="grey90"),
    axis.ticks = element_line(),
    text=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    legend.title=element_blank() ,#element_text(size=14, colour="black", face='plain'),
    plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm"))+
  theme(aspect.ratio=1)

ggsave(plot = pca_plot_ruv_type, file="pca_plot_ruv_type.svg", device=svg,
       bg="transparent", width=5, height=4, units=c("in"))