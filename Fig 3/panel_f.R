## Depends on Fig 3/panel_b.R and bulk/limma.R

# import y and meta from limma.R
y <- readRDS("./../bulk/dge_final/y.rds")
meta <- readRDS("./../bulk/dge_final/meta.rds")

# subset cpm count table
lcpm_subset_bulk <- cbind(
  edgeR::cpm(y)[,"M25"], #NSIP
  edgeR::cpm(y)[,"M29"],
  edgeR::cpm(y)[,"M33"],
  
  edgeR::cpm(y)[,"M3"], #CHP
 # edgeR::cpm(y)[,"M32"],
  edgeR::cpm(y)[,"M35"],
 # edgeR::cpm(y)[,"M38"],
  #edgeR::cpm(y)[,"M44"],
  
  edgeR::cpm(y)[,"M24"], #IPF
 # edgeR::cpm(y)[,"M34"],
  edgeR::cpm(y)[,"M37"]
  
)
#scale
lcpm_subset_bulk_scale <- t(scale(t(lcpm_subset_bulk)))

# import gene_order to match row/gene order of Nanostring heatmaps
gene_order_xfibrosis <- readRDS("./../nanostr/dge_final/across_fibrosis/gene_order_xfibrosis.rds")


# match returns indices in bulk topTags for top 50 genes from Nanostring
lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis <- lcpm_subset_bulk_scale[
  na.omit(match(gene_order_xfibrosis, rownames(lcpm_subset_bulk_scale))),
]


# complexheatmap
col_fun <- colorRamp2(c(-2,0,2), hcl_palette = 'Inferno')



chm.legend <- Legend(col_fun=col_fun,
                     border="black",
                     title="Z score", title_position="leftcenter-rot",
                     title_gp = gpar(fontsize=10,fontface='plain', fontfamily='sans'),
                     labels_gp = gpar(fontsize=10, fontface='plain', fontfamily='sans'),
                     legend_height=unit(1.8, "cm"),
                     legend_width=unit(5,"mm"))

svglite("legend.svg", width=1, height=2, bg="transparent")
draw(chm.legend)
dev.off()




chm_xfibrosis <- Heatmap(lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis, 
               col = col_fun,
               show_heatmap_legend = F,
               # heatmap_legend_param = list(
               #   #lgd = Legend(col_fun = col_fun, border='black')
               #   border="black",
               #   title="Z score", title_position="topcenter",
               #   title_gp = gpar(fontsize=9,fontface='plain', fontfamily='sans'),
               #   labels_gp = gpar(fontsize=8, fontface='plain', fontfamily='sans'),
               #   grid_height=unit(7, "mm"),
               #   grid_width=unit(3,"mm")
               # ),
               width= unit(dim(lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis)[2]*1.5, "mm"),
               height = unit(dim(lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis)[1]*1.5, "mm"),
               #right_annotation=ha,
               top_annotation = HeatmapAnnotation(
                 foo= anno_block(gp = gpar(lty=0, fill="transparent"),
                                 labels = c("NSIP","CHP","IPF"),
                                 labels_gp = gpar(col="black", fontsize=7, fontfamily='sans',fontface='bold'),
                                 labels_rot=45, labels_just = "center", labels_offset = unit(3,"mm"))
               ),
               
               # cluster_rows = dxt, row_dend_gp = gpar(lwd=1), row_split = 2, row_title=NULL,
               column_split = rep(LETTERS[1:3], times=c(3,2,2)) , column_title= NULL, column_gap = unit(0.2,"mm"),
               row_split = rep(LETTERS[1:2], times=c(38,12)), row_title= NULL, row_gap = unit(0.5,"mm"),
               border_gp =  gpar(col="black", lwd=0.2),
               show_column_names = F,
               row_names_gp = gpar(fontfamily = 'sans', fontface = 'italic', fontsize = 5),
               row_dend_width = unit(0.2, "cm"),
               cluster_columns=F, cluster_rows=F)

svg("chm_xfibrosis.svg", width=1.5, height=4, bg = "transparent")
draw(chm_xfibrosis, background="transparent")
dev.off()

