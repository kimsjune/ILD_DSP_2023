library(limma)
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(egg)
library(grid)
library(stringr)
library(ggrepel)
library(ggpp)
library(stringr)
library(scales)
library(cowplot)
library(ggpubr)
library(volcano3D)
library(ggfortify)
library(ggdendro)
library(dendextend)
library(scico)
library(standR)
library(SpatialExperiment)
library(pals)
library(ggpattern)
library(tidyverse)
library(complexheatmap)
library(circlize)

# contrasts of interest

fit <- readRDS(file= "./../fit.rds")
design <- readRDS(file= "./../design.rds")
dge_all <- readRDS(file="./../dge_all.rds")


con <- makeContrasts(
  fibrosis_CHP_v_central_NSIP = fibrosis_CHP - central_NSIP,
  fibrosis_IPF_v_fibroblast_IPF = fibrosis_IPF - fibroblast_IPF,
  fibrosis_IPF_v_central_NSIP = fibrosis_IPF - central_NSIP,
  fibrosis_CHP_v_fibroblast_IPF = fibrosis_CHP - fibroblast_IPF,
  fibrosis_CHP_v_fibrosis_IPF = fibrosis_CHP - fibrosis_IPF,
  fibroblast_IPF_v_central_NSIP =  fibroblast_IPF - central_NSIP,


  levels=design)


### following limma-voom pipeline 
set.seed(0)
fit_contrast <- contrasts.fit(fit, contrasts = con)
efit <- eBayes(fit_contrast, robust = TRUE)

# save efit to compare with bulk RNA-seq data
saveRDS(efit,"efit_across_fibrosis.rds")
topGenes <- topTable(efit, coef=c(4,6), n=50,   p.value=0.05, adjust.method="BH", lfc=1)


## RUV4-normalized log CPM whose columns are grouped by annotation

lcpm_subset <- cbind(
  assay(spe_ruv_subset,2)[, colData(spe_ruv_subset)$anno_type=="central_NSIP"],
  assay(spe_ruv_subset,2)[, colData(spe_ruv_subset)$anno_type=="fibrosis_CHP"],
  assay(spe_ruv_subset,2)[, colData(spe_ruv_subset)$anno_type=="fibrosis_IPF"],
  assay(spe_ruv_subset,2)[, colData(spe_ruv_subset)$anno_type=="fibroblast_IPF"]
  )

# scale
lcpm_subset_scale <- t(scale(t(lcpm_subset)))

# select genes
lcpm_subset_scale_topGenes <- lcpm_subset_scale[rownames(topGenes),]



# for complexheatmap
col_fun <- colorRamp2(c(-2,0,2), hcl_palette = 'Inferno')



chm.legend <- Legend(col_fun=col_fun,
                     border="black",
                     title="Z score", title_position="leftcenter-rot",
                     title_gp = gpar(fontsize=10,fontface='plain', fontfamily='sans'),
                     labels_gp = gpar(fontsize=10, fontface='plain', fontfamily='sans'),
                     legend_height=unit(1.8, "cm"),
                     legend_width=unit(5,"mm"))

svg("legend.svg", width=1, height=2, bg="transparent")
draw(chm.legend)
dev.off()



## row dendrogram
dendro <- as.dendrogram(hclust(d=dist(x=lcpm_subset_scale_topGenes)))


# maybe use dendextend instead of ggplot?
dxt <- dendro %>% 
  dendextend::set("branches_k_color", value=c("#c85200","#366785"), k=2) %>%
  dendextend::set("labels"," ")  %>%
  dendextend::set("branches_lwd", 0.5) 



chm <- Heatmap(lcpm_subset_scale_topGenes, 
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
        width= unit(dim(lcpm_subset_scale_topGenes)[2]*1.5, "mm"),
        height = unit(dim(lcpm_subset_scale_topGenes)[1]*1.5, "mm"),
        #right_annotation=ha,
        top_annotation = HeatmapAnnotation(
          foo= anno_block(gp = gpar(lty=0, fill="transparent"), 
                          labels = c("Central F. NSIP", "Fibrosis CHP", "Fibrosis IPF", "F.foci IPF"),
                          labels_gp = gpar(col="black", fontsize=7, fontfamily='sans', fontface='bold'),
                          labels_rot=20, labels_just = "center", labels_offset = unit(4,"mm"))
        ),
        
        cluster_rows = dxt, row_dend_gp = gpar(lwd=0.5), row_split = 2, row_title=NULL,
        column_split = rep(LETTERS[1:4], times=rep(9,4)) , column_title= NULL, column_gap = unit(0.2,"mm"),
        border_gp =  gpar(col="black", lwd=0.2),
        show_column_names = F,
        row_names_gp = gpar(fontfamily = 'sans', fontface = 'italic', fontsize = 5),
        row_dend_width = unit(0.2, "cm"), 
         cluster_columns=F)

svg("chm.svg", width=4, height=4, bg = "transparent")
draw(chm, background="transparent")
dev.off()


# extract row dendrogram to order genes for bulk tissue

gene_order <- c(rownames(lcpm_subset_scale_topGenes[ row_order(chm)[[1]],]),
                rownames(lcpm_subset_scale_topGenes[ row_order(chm)[[2]],])
                
  
)

saveRDS(gene_order, "gene_order_xfibrosis.rds")