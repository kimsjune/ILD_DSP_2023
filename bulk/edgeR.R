library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(egg)
library(grid)
library(stringr)
library(ggrepel)
library(stringr)
library(scales)

# import read count tables
samples <- c("M3","M24","M25","M26","M28","M29",
             "M30","M31","M32","M33",
             "M34","M35","M36","M37","M38",
             "M40","M41","M42","M43","M44")

M3 <- read.delim("./table/M3.csv", header=F, row.names=1)
M24 <- read.delim("./table/M24.csv", header=F, row.names=1)
M25 <- read.delim("./table/M25.csv", header=F, row.names=1)
M26 <- read.delim("./table/M26.csv", header=F, row.names=1)
M28 <- read.delim("./table/M28.csv", header=F, row.names=1)
M29 <- read.delim("./table/M29.csv", header=F, row.names=1)
M30 <- read.delim("./table/M30.csv", header=F, row.names=1)
M31 <- read.delim("./table/M31.csv", header=F, row.names=1)
M32 <- read.delim("./table/M32.csv", header=F, row.names=1)
M33 <- read.delim("./table/M33.csv", header=F, row.names=1)
M34 <- read.delim("./table/M34.csv", header=F, row.names=1)
M35 <- read.delim("./table/M35.csv", header=F, row.names=1)
M36 <- read.delim("./table/M36.csv", header=F, row.names=1)
M37 <- read.delim("./table/M37.csv", header=F, row.names=1)
M38 <- read.delim("./table/M38.csv", header=F, row.names=1)
M40 <- read.delim("./table/M40.csv", header=F, row.names=1)
M41 <- read.delim("./table/M41.csv", header=F, row.names=1)
M42 <- read.delim("./table/M42.csv", header=F, row.names=1)
M43 <- read.delim("./table/M43.csv", header=F, row.names=1)
M44 <- read.delim("./table/M44.csv", header=F, row.names=1)


df <- data.frame(M3, M24, M25, M26, M28, M29, M30, M31, M32, M33,
                 M34, M35, M36, M37, M38, M40, M41, M42, M43, M44)

df <- df[1:(nrow(df)-5),]
colnames(df) <- c("CHP","IPF","NSIP","UNC","COP",
                        "NSIP","IPAF","IPAF","CHP","NSIP",
                  "IPF","CHP","UNC","IPF","CHP","UNC",
                  "UNC","IPAF","UNC","CHP")

# edgeR and glm model
group <- factor(colnames(df))

y <- DGEList(df, group=group)


design <- model.matrix(~0+group)
colnames(design) <- levels(group)

keep <- filterByExpr(y, design)
table(keep)

y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)

y <- estimateDisp(y,design,robust=T)

fit <- glmFit(y, design, robust=T)

con <- makeContrasts(
  IPF_v_CHP = IPF - CHP,
  IPF_v_NSIP = IPF - NSIP,
  CHP_v_NSIP = CHP - NSIP,
  IPF_v_IPAF = IPF - IPAF,
  
  levels=design
)

lrt <- glmLRT(fit, contrast=con)
topTags(lrt)
summary(decideTests(lrt))
