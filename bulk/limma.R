## Using limma-voom to see if any genes that were identified as differential in spatial data are also differetially expressed in bulk tissues
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(egg)
library(grid)
library(stringr)
library(ggrepel)
library(stringr)
library(scales)
library(limma)

packages <- c("edgeR","limma","standR","SpatialExperiment",
              "egg","grid","cow",
              "scales","stringr",
              "ggplot2","ggpubr","ggrepel","ggpp",
              "ggdendro","dendextend",
              "tidyverse")
lapply(packages, require, character.only = TRUE)

### import count tables

samples <- c("M3","M24","M25","M29","M32","M33",
             "M34","M35","M37","M38","M44")

M3 <- read.delim("./counts/M3.csv", header=F, row.names=1)
M24 <- read.delim("./counts/M24.csv", header=F, row.names=1)
M25 <- read.delim("./counts/M25.csv", header=F, row.names=1)
M29 <- read.delim("./counts/M29.csv", header=F, row.names=1)
M32 <- read.delim("./counts/M32.csv", header=F, row.names=1)
M33 <- read.delim("./counts/M33.csv", header=F, row.names=1)

M34 <- read.delim("./counts/M34.csv", header=F, row.names=1)
M35 <- read.delim("./counts/M35.csv", header=F, row.names=1)
M37 <- read.delim("./counts/M37.csv", header=F, row.names=1)
M38 <- read.delim("./counts/M38.csv", header=F, row.names=1)
M44 <- read.delim("./counts/M44.csv", header=F, row.names=1)


### pool then remove last few non gene rows
df <- data.frame(M3, M24, M25, M29, M32, M33,
                 M34, M35, M37, M38, M44)
df <- df[1:(nrow(df)-5),]
# unique column names
colnames(df) <- samples


### create sample metadata table
# patid, type, sequencing batch

meta <- data.frame(
  patid = samples,
  type = c("CHP","IPF","NSIP","NSIP","CHP","NSIP",
           "IPF","CHP","IPF","CHP","CHP"),
  batch = c(rep(1,6),rep(2,5)),
  replicate = c("1","1","1","2","2","3",
                "2","3","3","4","5")
)

saveRDS(meta,"meta.rds")
group <- factor(meta$type)
batch <- factor(meta$batch)

# limma to find DGEs

y <- DGEList(df, group=group)

y <- calcNormFactors((y))



# model matrix, sequencing batch # as confounder

design <- model.matrix(~0+group+batch)
colnames(design) <- gsub("group","",colnames(design))

keep <- filterByExpr(y, design)


y <- y[keep,,keep.lib.sizes=FALSE]
y <- estimateDisp(y, design, robust=T)

plotMDS(y, col = as.numeric(meta$type))




# variance mean

v <- voom(y, design, plot = T)

# create contrast matrix

con <- makeContrasts(
  IPF_v_CHP = IPF - CHP,
  IPF_v_NSIP = IPF - NSIP,
  CHP_v_NSIP = CHP - NSIP,

  
  levels=design
)


# fit linear model

fit <- glmFit(y, design, robust=T)
fit.limma <- lmFit(v, design, robust=T)
head(coef(fit))

# # LR test
# 
# lrt <- glmLRT(fit, contrast=con)
# 
# 
# topTags(lrt)



fit.limma <- contrasts.fit(fit.limma, contrasts=con)
efit <- eBayes(fit.limma)

topTable(efit, coef=c(1:3))
summary_efit <-summary(decideTests(efit))

saveRDS(efit,"fit.rds")
saveRDS(design,"design.rds")
saveRDS(y,"y.rds")


