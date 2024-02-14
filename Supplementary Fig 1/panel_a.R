## ROI nuclei count vs. Library size plot, grouped by condition
# depends on standR/import.R 
# customized plotROIQC() function from standR package

plotROIQC_custom <- function(spe_object,
                      
                      x_axis = "AOINucleiCount",
                      y_axis = "lib_size",
                      x_lab = "AOINucleiCount",
                      y_lab = "Library size",
                      x_threshold = NULL,
                      y_threshold = NULL,
                      regression_col = "purple",
                      hist_col = "black", hist_fill = "white", bin_num = 50,
                      threshold_col = "red", threshold_linetype = "dashed",
                      layout_ncol = 2, layout_nrow = 2,
                      leyout_height = c(0.8, 2.5), layout_width = c(2.5, 0.8),
                      ...) {
  stopifnot(x_axis %in% colnames(colData(spe_object)))
  stopifnot(y_axis %in% colnames(colData(spe_object)))
  
  aesmap <- rlang::enquos(...)
  x_axis <- rlang::sym(x_axis)
  y_axis <- rlang::sym(y_axis)
  
  # plot dot plot
  p1 <- colData(spe_object) |>
    as.data.frame(optional = TRUE) |>
    ggplot(aes(!!x_axis, !!y_axis, !!!aesmap)) +
    geom_point(aes(fill=factor(type)),colour="black", pch=21, size=2.5)+
    
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

    # geom_smooth(method = "loess", se = FALSE, col = regression_col) +
    #theme_test() +
    xlab(x_lab) +
    ylab(y_lab) +
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
    guides(fill=guide_legend(nrow=2, byrow=T))+
    xlab("ROI nuclei count")
  
  # plot distribution of y axis
  p2 <- colData(spe_object) |>
    as.data.frame(optional = TRUE) |>
    ggplot(aes(!!y_axis)) +
    geom_histogram(col = hist_col, fill = hist_fill, bins = bin_num) +
    #theme_test() +
    coord_flip() +
    theme_bw()+
    theme(text=element_text(size=14),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size=14, color="black"),
          panel.grid.minor=element_blank()                           
    )+
    ylab("Frequency") 
  
  p_blank <- ggplot() +
    theme_void()
  
  # plot distribution of x axis
  p3 <- colData(spe_object) |>
    as.data.frame(optional = TRUE) |>
    ggplot(aes(!!x_axis)) +
    geom_histogram(col = hist_col, fill = hist_fill, bins = bin_num) +
    #theme_test() +
    theme_bw()+
    theme(text=element_text(size=14),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=14, color="black"),
          panel.grid.minor=element_blank()  )+
    ylab("Frequency") 
  
  # plot threshold
  if (!is.null(x_threshold)) {
    stopifnot(is.numeric(x_threshold))
    p1 <- p1 + geom_vline(xintercept = x_threshold, col = threshold_col, linetype = threshold_linetype)
    p3 <- p3 + geom_vline(xintercept = x_threshold, col = threshold_col, linetype = threshold_linetype)
  }
  
  if (!is.null(y_threshold)) {
    stopifnot(is.numeric(y_threshold))
    p1 <- p1 + geom_hline(yintercept = y_threshold, col = threshold_col, linetype = threshold_linetype)
    p2 <- p2 + geom_vline(xintercept = y_threshold, col = threshold_col, linetype = threshold_linetype)
  }
  
  # p3 + p_blank + p1 + p2 + patchwork::plot_layout(layout_ncol, layout_nrow,
  #                                                 widths = layout_width,
  #                                                 heights = leyout_height, guides = "collect")
  p1
  
}

ROI_plot <- plotROIQC_custom(speGQC, color=type,
                             y_threshold=NULL) 
ROI_plot_resized <- set_panel_size(ROI_plot,
  width=unit(3,"in"), height=unit(3,"in"), margin = margin(0.2, 0.2, 0.2, 0.2)
)
ggsave(plot= ROI_plot_resized, file="ROI_plot.svg", bg="transparent", device=svg,
       width=5, height=5, units=c("in"))

       ## ROI nuclei vs. Library size plot, grouped by annotation
       plotROIQC_custom2 <- function(spe_object,
                             
                             x_axis = "AOINucleiCount",
                             y_axis = "lib_size",
                             x_lab = "AOINucleiCount",
                             y_lab = "Library size",
                             x_threshold = NULL,
                             y_threshold = NULL,
                             regression_col = "purple",
                             hist_col = "black", hist_fill = "white", bin_num = 50,
                             threshold_col = "red", threshold_linetype = "dashed",
                             layout_ncol = 2, layout_nrow = 2,
                             leyout_height = c(0.8, 2.5), layout_width = c(2.5, 0.8),
                             ...) {
  stopifnot(x_axis %in% colnames(colData(spe_object)))
  stopifnot(y_axis %in% colnames(colData(spe_object)))
  
  aesmap <- rlang::enquos(...)
  x_axis <- rlang::sym(x_axis)
  y_axis <- rlang::sym(y_axis)
  
  # plot dot plot
  p1 <- colData(spe_object) |>
    as.data.frame(optional = TRUE) |>
    ggplot(aes(!!x_axis, !!y_axis, !!!aesmap)) +
    geom_point(aes(fill=factor(anno)),colour="black", pch=21, size=2.5)+
    
    #geom_point(size=3) +
    scale_fill_manual("Annotation",
      guide="legend",
                      values=c(
                      "airway"="#c10534",
                        "central"="#938dd2",
                        "fibroblast"="#cac27e",
                        "fibrosis"="#a0522d",
                        "granuloma"="#7b92a8",
                        "inflammatory"="#2d6d66",
                        "lymphoid"="#9c8847",
                        "neutral"="#bfa19c",
                        "peripheral"="#ffd200",
                        "pleura"="#d9e6eb"
                      ),
                      breaks=c("airway","central","fibroblast","fibrosis","granuloma",
                               "inflammatory","lymphoid","neutral","peripheral",
                               "pleura"),
                      labels=c("airway","central","FF","fibrosis","granuloma",
                               "inflammatory","lymphoid","uninvolved","peripheral",
                               "pleura"))+
    
   # geom_smooth(method = "loess", se = FALSE, col = regression_col) +
    #theme_test() +
    xlab(x_lab) +
    ylab(y_lab) +
    theme_bw() +
    theme(legend.text=element_text(size=14, colour="black"),
          legend.title=element_text(size=14, colour="black"),
          text=element_text(size=14, color="black"),
          axis.text=element_text(size=14, color="black"),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin=unit(c(0.2, 0.6, 0.2, 0.2), "cm"),
          aspect.ratio=1,
          plot.background = element_blank())+
   # guides(fill=guide_legend(nrow=3, byrow=T))+
    xlab("ROI nuclei count")
  
  # plot distribution of y axis
  p2 <- colData(spe_object) |>
    as.data.frame(optional = TRUE) |>
    ggplot(aes(!!y_axis)) +
    geom_histogram(col = hist_col, fill = hist_fill, bins = bin_num) +
    #theme_test() +
    coord_flip() +
    theme_bw()+
    theme(text=element_text(size=14),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size=14, color="black"),
          panel.grid.minor=element_blank()                           
    )+
    ylab("Frequency") 
  
  p_blank <- ggplot() +
    theme_void()
  
  # plot distribution of x axis
  p3 <- colData(spe_object) |>
    as.data.frame(optional = TRUE) |>
    ggplot(aes(!!x_axis)) +
    geom_histogram(col = hist_col, fill = hist_fill, bins = bin_num) +
    #theme_test() +
    theme_bw()+
    theme(text=element_text(size=14),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=14, color="black"),
          panel.grid.minor=element_blank()  )+
    ylab("Frequency") 
  
  # plot threshold
  if (!is.null(x_threshold)) {
    stopifnot(is.numeric(x_threshold))
    p1 <- p1 + geom_vline(xintercept = x_threshold, col = threshold_col, linetype = threshold_linetype)
    p3 <- p3 + geom_vline(xintercept = x_threshold, col = threshold_col, linetype = threshold_linetype)
  }
  
  if (!is.null(y_threshold)) {
    stopifnot(is.numeric(y_threshold))
    p1 <- p1 + geom_hline(yintercept = y_threshold, col = threshold_col, linetype = threshold_linetype)
    p2 <- p2 + geom_vline(xintercept = y_threshold, col = threshold_col, linetype = threshold_linetype)
  }
  
  # p3 + p_blank + p1 + p2 + patchwork::plot_layout(layout_ncol, layout_nrow,
  #                                                 widths = layout_width,
  #                                                 heights = leyout_height, guides = "collect")
  # 
  p1
}

ROI_plot_anno <- plotROIQC_custom2(speGQC, color=anno,
                             y_threshold=NULL) 

ROI_plot_anno_resized <- set_panel_size(
  ROI_plot_anno, width=unit(3,"in"), height=unit(3,"in"),
  margin = margin(0.2, 0.2, 0.2, 0.2)
)
ggsave(plot= ROI_plot_anno_resized, file="ROI_plot_anno.svg", bg="transparent", device=svg,
       width=7, height=5, units=c("in"))
#########################