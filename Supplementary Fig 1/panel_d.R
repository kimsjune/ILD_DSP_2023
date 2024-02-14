## PCA with upperquartile norm

spe_q3 <- geomxNorm(speFinal, method="upperquartile")
spe_q3 <- scater::runPCA(spe_q3)
pca_results_q3 <- reducedDim(spe_q3, "PCA")




## by patient id
pca_plot_q3 <- drawPCA(spe_q3, precomputed = pca_results_q3) +
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

ggsave(plot = pca_plot_q3, file="pca_plot_q3.svg", device=svg,
       bg="transparent", width=5, height=4, units=c("in"))

## by type 
pca_plot_q3_type <- drawPCA(spe_q3, precomputed = pca_results_q3) +
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

ggsave(plot = pca_plot_q3_type, file="pca_plot_q3_type.svg", device=svg,
       bg="transparent", width=5, height=4, units=c("in"))





# make type into a factor
colData(spe_q3)$type <- factor(colData(spe_q3)$type, 
                               levels=c("IPF","NSIP","CHP","UNC","NOR"),
                               ordered=T)
# returns colData where rows are sorted by type from IPF to NOR
colData(spe_q3)[do.call(order, colData(spe_q3)["type"]),]

## Relative log2 expression

RLExpr_sizeNorm <- plotRLExpr(
  spe_q3[,rownames(colData(spe_q3)[do.call(order, colData(spe_q3)["type"]),])], 
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

ggsave(plot= RLExpr_sizeNorm, file="RLExpr_sizeNormQ3.svg", bg="transparent", device=svg,
       width=5, height=4, units=c("in"))