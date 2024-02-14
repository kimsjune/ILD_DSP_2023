## Relative log2 expression

RLExpr_cpm <- plotRLExpr(
  # rownames of colData(...) is used to sort the columns of speFinal 
  # no amount of using 'breaks' to indicate order of types works
  speFinal[,rownames(colData(speFinal)[do.call(order, colData(speFinal)["type"]),])],
                         assay=2, color=type, ordannots="type")+
  geom_boxplot()+
  geom_point(colour="black",size=1)+
  #scale_x_discrete(breaks=c("IPF","NSIP","CHP","UNC","NOR"))+
  scale_colour_manual(
    values=c(
      "CHP"="#1a476f",
      "IPF"= "#90353b",
      "NOR"="#55752f",
      "NSIP"="#e37e00",
      "UNC" ="#6e8e84"
    ))+
  xlab("ROIs")+

  theme(axis.ticks=element_line(),
        legend.text=element_text(size=14, colour="black"),
        legend.title=element_blank(),
        text=element_text(size=14, color="black"),
        axis.text=element_text(size=14, color="black"),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm"))+
  theme(aspect.ratio=1)

ggsave(plot= RLExpr_cpm, file="RLExpr_cpm.svg", bg="transparent", device=svg,
       width=5, height=4, units=c("in"))


## PCA

speFinal <- scater::runPCA(speFinal)
pca_results_cpm <- reducedDim(speFinal, "PCA")

# by patient id
pca_plot_cpm <- drawPCA(speFinal, precomputed = pca_results_cpm) +
  # scale_color_manual(values=c(stata_pal("s2color")(16)))
  # scale_color_tableau(
  #   palette="Classic 20",
  #   type=c("ordered-diverging"),
  #   direction=1
  geom_point(aes(fill=factor(patid)),colour="black", pch=21, size=2.5)+
  scale_fill_manual(guide="legend",
                     
                      values=c(unname(glasbey())),
                      breaks=c(mixedsort(unique(speFinal$patid))))+
  # scale__discrete(breaks=c(mixedsort(unique(speFinal$patid))))+
  theme(#panel.background = element_rect(fill="grey90"),
        axis.ticks = element_line(),
        text=element_text(size=14, color="black"),
        axis.text=element_text(size=14, color="black"),
        legend.title=element_blank() ,#element_text(size=14, colour="black", face='plain'),
        plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm"))+
  theme(aspect.ratio=1)
ggsave(plot = pca_plot_cpm, file="pca_plot_cpm.svg", device=svg,
       bg="transparent", width=5, height=4, units=c("in"))


# by disease

pca_plot_cpm_type <- drawPCA(speFinal, precomputed = pca_results_cpm) +
  scale_fill_manual(
    values=c(
      "CHP"="#1a476f",
      "IPF"= "#90353b",
      "NOR"="#55752f",
      "NSIP"="#e37e00",
      "UNC" ="#6e8e84"
    ),
    breaks=c("IPF","NSIP","CHP","UNC","NOR"))+
  
  geom_point(aes(fill=factor(type)), colour="black", pch=21, size=2.5)+
  # scale_color_manual(values=c(stata_pal("s2color")(16)))
  # scale_color_tableau(
  #   palette="Classic 20",
  #   type=c("ordered-diverging"),
  #   direction=1
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

ggsave(plot = pca_plot_cpm_type, file="pca_plot_cpm_type.svg", device=svg,
       bg="transparent", width=5, height=4, units=c("in"))