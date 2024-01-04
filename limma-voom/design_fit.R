### dge with ruv

colData(spe_ruv)[,seq(ncol(colData(spe_ruv))-1, ncol(colData(spe_ruv)))] |>
  head()


library(edgeR)
library(limma)

dge <- SE2DGEList(spe_ruv)

design <- model.matrix(~0 + anno_type + ruv_W1 + ruv_W2 + ruv_W3 + ruv_W4,
                       data = colData(spe_ruv))
colnames(design) <- gsub("anno_type","", colnames(design))
colnames(design)



keep <- filterByExpr(dge, design)
table(keep)
dge_all <- dge[keep, , keep.lib.sizes=F]


### BCV check
dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)

# plot 
plotBCV(dge_all, legend.position = "topleft", ylim = c(0, 1.3))
bcv_df <- data.frame(
  'BCV' = sqrt(dge_all$tagwise.dispersion),
  'AveLogCPM' = dge_all$AveLogCPM,
  'gene_id' = rownames(dge_all)
)

highbcv <- bcv_df$BCV > 0.8
highbcv_df <- bcv_df[highbcv, ]
points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "red")
text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4)


# differential expression
# double voom and dupliateCorrelation
# first
v <- voom(dge_all, design, plot = TRUE) 
corfit <- duplicateCorrelation(v, design, block=colData(spe_ruv)$patid)

# second
v <- voom(dge_all, design, plot=T, block=colData(spe_ruv)$patid, correlation=corfit$consensus)
corfit2 <- duplicateCorrelation(v, design, block=colData(spe_ruv)$patid)

fit <- lmFit(v, design, block=colData(spe_ruv)$patid, correlation=corfit2$consensus)

saveRDS(fit, "fit.rds")
saveRDS(design, "design.rds")
saveRDS(dge_all, "dge_all.rds")


