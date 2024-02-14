# % sequencing saturation vs. alignment plot

alignment_saturation <- ggplot(data=as.data.frame(colData(speGQC)),
                          aes(SequencingSaturation,AlignedReads*100/RawReads))+
  geom_point(aes(fill=factor(type)),colour="black", pch=21, size=2.5)+
  geom_vline(xintercept = 90, col = "red", linetype = "dashed")+
  geom_hline(yintercept=90, col="red", linetype="dashed")+
  
  #geom_point(size=3) +
  scale_fill_manual("Condition",
                    guide="legend",
                    values=c(
                      "CHP"="#1a476f",
                      "IPF"= "#90353b",
                      "NOR"="#55752f",
                      "NSIP"="#e37e00",
                      "UNC" ="#6e8e84"
                    ),
                    breaks=c("IPF","NSIP","CHP","UNC","NOR"))+
  scale_x_continuous(limits=c(0,100))+
  scale_y_continuous(limits=c(0,100))+
  
  #geom_smooth(method = "loess", se = FALSE, col = regression_col) +
  #theme_test() +
  xlab("% Sequencing Saturation") +
  ylab("% Alignment") +
  theme_bw() +
  theme(legend.text=element_text(size=14, colour="black"),
        legend.title=element_text(size=14, colour="black"),
        legend.background=element_rect(fill="transparent", colour=NA),
        legend.box.background=element_rect(fill="transparent", colour=NA),
        text=element_text(size=14, color="black"),
        axis.text=element_text(size=14, color="black"),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        aspect.ratio=1,
        legend.position = "top",
        plot.background = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=T))
  #xlab("ROI nuclei count")

alignment_saturation_resized <- set_panel_size(alignment_saturation,
                                   width=unit(3,"in"), height=unit(3,"in"), margin = margin(0.2, 0.2, 0.2, 0.2)
)

ggsave(plot=alignment_saturation_resized, file="alignment_saturation.svg",
       device=svg, bg="transparent",
       width=5, height=5, units=c("in"))



### remove ROIs with low % alignement OR sequencing saturation

qc <- colData(speGQC)$AlignedReads/colData(speGQC)$RawReads >=0.9 & colData(speGQC)$SequencingSaturation >=90
table(qc)
speFinal <- speGQC[,qc]
