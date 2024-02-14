## alluvial plot

alluvial_plot <- colData(speFinal) |>
  as.data.frame(optional=TRUE) |>
  dplyr::select(c("type","anno")) |>
  dplyr::mutate(anno = gsub("neutral","uninvolved", anno))|>
  dplyr::mutate(anno = gsub("fibroblast","fibroblastic foci", anno))|>
  ggalluvial::to_lodes_form() |>
  ggplot(aes(
    x = x, stratum = factor(stratum, 
                            levels=c("IPF","NSIP","CHP","UNC","NOR",
                                     "airway","central","fibroblastic foci","fibrosis",
                                     "granuloma","inflammatory","lymphoid",
                                     "uninvolved","peripheral","pleura")),
    alluvium = alluvium, fill = stratum, label = stratum, alpha=0.8
  )) +
  ggalluvial::geom_flow(
    stat = "alluvium", lode.guidance = "frontback",
    color = "NA"
  ) +
  theme_bw()+
  ggalluvial::geom_stratum(alpha=0.7)+
  geom_text(stat="stratum", color="black", size=5, alpha=1,
            aes(
              label = after_stat(stratum),
              hjust = ifelse(x == "type", 1, 0),
              x = as.numeric(factor(x)) + .075 * ifelse(x == "anno", 3, -3)
            ))+
  scale_fill_manual(values=c(stata_pal("s2color")(15)))+
  scale_y_discrete(breaks=c("IPF","NSIP","CHP","UNC","NOR"
  ))+
  scale_x_discrete(labels=c("Condition","Annotation"),
                   expand=expansion(add=c(0.5,0.75)))+
  scale_y_continuous(name="Frequency")+
  theme(legend.position="none",
        text=element_text(size=14, color="black"),
        axis.text=element_text(size=14, color="black"),
        axis.title.x=element_blank())


ggsave(plot=alluvial_plot,
       filename="alluvial_plot.svg", bg="transparent", device=svg,
       width=6.5, height=4, units=c("in"))