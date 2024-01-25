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


### using limma voom fit model and design
# see standR_RUV_limma_voom.R

fit <- readRDS(file= "./../fit.rds")
design <- readRDS(file= "./../design.rds")
dge_all <- readRDS(file="./../dge_all.rds")

con <- makeContrasts(
  fibroblast_v_fibrosis = fibroblast_IPF -fibrosis_IPF,
  fibrosis_v_neutral = fibrosis_IPF - neutral_IPF,
  fibroblast_v_neutral = fibroblast_IPF - neutral_IPF,
  #fibrosis_v_lymphoid = fibrosis_IPF - lymphoid_IPF,
  #lymphoid_v_neutral = lymphoid_IPF - neutral_IPF,
  levels=colnames(design))

