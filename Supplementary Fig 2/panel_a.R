####
### reverse deconvolution with SpatialDecon example
####

## 1:fibroblast_v_fibrosis subset of spe data
subset_IPF_fibroblast_v_fibrosis <- c(
  colnames(spe_keepNeg_ruv)[colData(spe_keepNeg_ruv)$anno_type == "neutral_IPF"],
  colnames(spe_keepNeg_ruv)[colData(spe_keepNeg_ruv)$anno_type == "fibroblast_IPF"]
)

# subset the res object from running spatialdecon
rdecon_IPF_fibroblast_v_fibrosis <- reverseDecon(
  norm = as.matrix(assay(spe_keepNeg_ruv,"logcounts")[,subset_IPF_fibroblast_v_fibrosis]),
  beta = res$beta[,subset_IPF_fibroblast_v_fibrosis]
  
)

# gene names to highlight in the biplot
showgenes_IPF_fibroblast_v_fibrosis = na.omit(readRDS("./../within_IPF/IPF_fibroblast_v_fibrosis_deLab.rds"))


# ggplot

df.IPF_fibroblast_v_fibrosis <- data.frame(name = names(rdecon_IPF_fibroblast_v_fibrosis$cors), x=rdecon_IPF_fibroblast_v_fibrosis$cors, y= rdecon_IPF_fibroblast_v_fibrosis$resid.sd)

# write table of x and y values
# write.table(df.IPF_fibroblast_v_fibrosis[showgenes_IPF_fibroblast_v_fibrosis,],
#             file = "IPF_fibroblast_v_fibrosis_xy.txt", sep = "\t", quote=F, row.names = F)
write.table(
  df.IPF_fibroblast_v_fibrosis%>%
  mutate(match = ifelse( name %in% showgenes_IPF_fibroblast_v_fibrosis,"DE","Not DE")) %>%
  arrange(match),
  "IPF_fibroblast_v_fibrosis_corr_v_sd.txt",
  sep = "\t", row.names =F, quote =F)

biplot.IPF_fibroblast_v_fibrosis <- ggplot(df.IPF_fibroblast_v_fibrosis, 
                                        aes(x=x, y=y)) + 
  theme_bw()+
  theme(
    axis.ticks = element_line(colour="black"),
    panel.border = element_rect(colour="black"),
    text=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    plot.margin=unit(c(1,0.1,1,0),"cm"),
    plot.background=element_rect(fill="transparent", colour=NA),
    panel.background=element_rect(fill="transparent", colour=NA),
    panel.grid = element_blank(),
    legend.background=element_rect(fill="transparent", colour=NA),
    legend.box.background=element_rect(fill="transparent", colour=NA),
    legend.key=element_rect(fill="transparent", colour=NA),
    legend.position="none")+
  geom_point(data = . %>%
               mutate(match = ifelse( name %in% showgenes_IPF_fibroblast_v_fibrosis,"DE","Not DE")),
             aes(color=match, alpha = match))+
  # geom_text_repel(data = . %>%
  #                   mutate(label = ifelse(name %in% showgenes_IPF_fibroblast_v_fibrosis, name, "")),
  #                 aes(label=label), direction = "both", xlim=c(0,1.5), box.padding = 0.5,
  #                 position = 
  #                   position_nudge_center(direction= "radial",
  #                                         center_x = 0.8,
  #                                         center_y= 0.1,
  #                                         x= 0.1, y= 0.05,
  #                                         ),
  #                 max.overlaps = 15)+
  scale_colour_manual(values = c("DE"="red4", "Not DE" = "grey50"))+
  scale_alpha_manual(values = c("DE"=1,"Not DE"=0.1))+
  # xlab(str_wrap("Correlation between predicted and observed expression",40))+
   ylab(str_wrap("SD of residuals from predicted expression",40))+
  xlab("")+
 # ylab("")+
  xlim(0,1)+
  ylim(0,0.25)+
  theme(aspect.ratio=1, rect=element_rect(fill="transparent"))


ggsave(plot=biplot.IPF_fibroblast_v_fibrosis,
       filename="biplot.IPF_fibroblast_v_fibrosis.svg", bg="transparent",
       width=4, height=4, units=c("in"), device=svg)