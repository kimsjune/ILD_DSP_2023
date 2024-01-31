## Output if genes are differentially expressed in bulk samples
fit <- readRDS("./../bulk/dge_final/fit.rds")
topTable(fit, coef=c(1:3), n=Inf)[rownames(lcpm_subset_bulk_scale_with_topGenesPathDistinctTypes),]


topTable(fit, coef=c(1:3), n=Inf)[rownames(lcpm_subset_bulk_scale_with_topGenesAcrossFibrosis),]