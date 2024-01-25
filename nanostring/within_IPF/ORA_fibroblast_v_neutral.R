# neutral refers to uninvolved

ego_fibroblast_v_neutral_up <-enrichGO(rownames(topT.fibroblast_v_neutral[topT.fibroblast_v_neutral$logFC > 1 & topT.fibroblast_v_neutral$adj.P.Val < 0.05,]), OrgDb=org.Hs.eg.db, ont="BP", keyType = "SYMBOL")

dot_fibroblast_v_neutral_up <- dotplot(ego_fibroblast_v_neutral_up, showCategory=3)+
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

dot_fibroblast_v_neutral_up.noLegend <- dot_fibroblast_v_neutral_up  + theme(legend.position="none")
dot_fibroblast_v_neutral_up.resized <- set_panel_size(
  dot_fibroblast_v_neutral_up.noLegend,
  width=unit(2,"in"), height=unit(1.5,"in")
)
ggsave(plot=dot_fibroblast_v_neutral_up.resized, 
       filename="dot_fibroblast_v_neutral_up.resized.svg", bg="transparent", device=svg,
       width=5, height=4, units=c("in"))

dot_fibroblast_v_neutral_up.legend <- get_legend(dot_fibroblast_v_neutral_up)
dot_fibroblast_v_neutral_up.legend <- as_ggplot(dot_fibroblast_v_neutral_up.legend)
ggsave(plot=dot_fibroblast_v_neutral_up.legend,
       filename="dot_fibroblast_v_neutral_up.legend.svg", bg="transparent",
       device=svg,
       width=6, height=1.5, units=c("in"))


#####

ego_fibroblast_v_neutral_dn <-enrichGO(rownames(topT.fibroblast_v_neutral[topT.fibroblast_v_neutral$logFC < -1 & topT.fibroblast_v_neutral$adj.P.Val < 0.05,]), OrgDb=org.Hs.eg.db, ont="BP", keyType = "SYMBOL")

dot_fibroblast_v_neutral_dn <- dotplot(ego_fibroblast_v_neutral_dn, showCategory=3)+
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

dot_fibroblast_v_neutral_dn.noLegend <- dot_fibroblast_v_neutral_dn  + theme(legend.position="none")
dot_fibroblast_v_neutral_dn.resized <- set_panel_size(
  dot_fibroblast_v_neutral_dn.noLegend,
  width=unit(2,"in"), height=unit(1.5,"in")
)
ggsave(plot=dot_fibroblast_v_neutral_dn.resized, 
       filename="dot_fibroblast_v_neutral_dn.resized.svg", bg="transparent", device=svg,
       width=5, height=4, units=c("in"))

dot_fibroblast_v_neutral_dn.legend <- get_legend(dot_fibroblast_v_neutral_dn)
dot_fibroblast_v_neutral_dn.legend <- as_ggplot(dot_fibroblast_v_neutral_dn.legend)
ggsave(plot=dot_fibroblast_v_neutral_dn.legend,
       filename="dot_fibroblast_v_neutral_dn.legend.svg", bg="transparent",
       device=svg,
       width=6, height=1.5, units=c("in"))

