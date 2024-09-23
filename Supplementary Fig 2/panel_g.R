## cnetplot
# depends on Fig 2/panel_g.R

up_cnet <- cnetplot(layout="kk",
  list("ECM organization"=
         c(unlist(strsplit(IPF_ego_fibroblast_v_fibrosis_up$geneID[[1]],"/")),
           unlist(strsplit(IPF_ego_fibroblast_v_neutral_up$geneID[[1]],"/"))),
       "Complement activation"=
         unlist(strsplit(IPF_ego_fibrosis_v_neutral_up$geneID[[1]],"/")),
       "Regulation of peptidase activity"=
         c(unlist(strsplit(NSIP_ego_central_v_inflammatory_up$geneID[[3]],"/")),
           unlist(strsplit(NSIP_ego_peripheral_v_inflammatory_up$geneID[[1]],"/"))
           
         )
)) + 
  theme(legend.position="none",
        plot.margin=unit(c(1,1,1,1),"cm"),
        plot.background=element_rect(fill="transparent", colour=NA),
        panel.background=element_rect(fill="transparent", colour=NA),
        panel.border=element_blank(),
        text=element_text(size=14, face='plain', family='sans'))

up_cnet.resized <- set_panel_size(
  up_cnet,
  width=unit(9,"in"), height=unit(6,"in")
)

ggsave(plot=up_cnet.resized,
       filename="up_cnet.resized.svg", bg='transparent', device=svg,
       width=9, height=6, units=c("in"))
