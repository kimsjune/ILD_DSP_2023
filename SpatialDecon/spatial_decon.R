####
### create a new spe object keeping all negative probes
####
# countFile has to be modified to keep negative probe rows. Columns 2 and 3 are kept
countFile_keepNeg <- as.data.frame(bioprob[,c(2,3,13:191)])

# Negative probes must have unique names - stored in the fourth column

# If it's NOT negative probes, don't keep V2, keep V3
notNeg <- countFile_keepNeg[countFile_keepNeg$V3!="NegProbe-WTX",-which(names(countFile_keepNeg) %in% c("V2"))]
# If it IS negative probes, dont keep V3, keep V2
Neg <-  countFile_keepNeg[countFile_keepNeg$V3=="NegProbe-WTX",-which(names(countFile_keepNeg) %in% c("V3"))]


# colnames must match to rbind
names(notNeg)[names(notNeg)=='V3'] <- 'V2'
combine <- rbind(notNeg,Neg)

# Under TargetName (V3), pasted ProbeNames that has numbers after NegProbe-WTX so they are all unique.
bioprob_keepNeg <- read.table("./../../dsp_download/bioprobe_remove_GLI_keepNeg.txt",
                              header=F,
                              sep="\t")
featureAnnoFile_keepNeg <- as.data.frame( bioprob_keepNeg[,c(1:12)])


featureAnnoFile_keepNeg <- rbind(
  featureAnnoFile_keepNeg[featureAnnoFile_keepNeg$V3!=grepl("NegProbe-WTX",featureAnnoFile_keepNeg$V3),],
  featureAnnoFile_keepNeg[!featureAnnoFile_keepNeg$V3!=grepl("NegProbe-WTX",featureAnnoFile_keepNeg$V3),])
#featureAnnoFile[featureAnnoFile$V3==grepl("NegProbe-WTX",featureAnnoFile$V3),])
# for some reason, thes tables must be written out as a file then imported directly into readGeoMx
write.table(combine,file="export_countFile_keepNeg.txt",sep="\t",
            row.names=F, col.names=F)
write.table(featureAnnoFile_keepNeg,file="export_featureAnnoFile_keepNeg.txt",sep="\t",
            row.names=F, col.names=F)




####
### Create spe and RUV normalize like negative probe removed data
####


spe_keepNeg <- readGeoMx("./export_countFile_keepNeg.txt",
                         "./export_sampleAnnoFilev2.txt",
                        "./export_featureAnnoFile_keepNeg.txt",
                         rmNegProbe = FALSE
                        # NegProbeName="NegProbe-WTX_137"
)


# QC steps

qc2 <- colData(spe_keepNeg)$AlignedReads/colData(spe_keepNeg)$RawReads >=0.9 & colData(spe_keepNeg)$SequencingSaturation >=90

spe_keepNeg <- spe_keepNeg[,qc2]

spe_keepNeg <- findNCGs(spe_keepNeg, batch_name="patid", top_n=200)

spe_keepNeg_ruv <- geomxBatchCorrection(spe_keepNeg, factors = "anno_type", 
                                        NCGs = metadata(spe_keepNeg)$NCGs, k = 4)


####
### Spatial Decon
####

# calculate background
bg <- derive_GeoMx_background(norm = assay(spe_keepNeg_ruv, "logcounts"),
                              probepool = rep(1, nrow(assay(spe_keepNeg_ruv, "logcounts"))),
                              negnames = paste0("NegProbe-WTX_",c(1:139)))


# download adult lung cell profile
lung <- download_profile_matrix(species = "Human", age_group = "Adult", matrixname = "Lung_HCA")

# too many cell types to plot individually, so group them into 10 bins
matching = list()
matching$myeloid = c("monocyte","MARCOneg.macrophage", "MARCOpos.macrophage")
matching$T.NK = c("CD4+.T.cell","CD8+.cytotoxic.T.cell", "regulatory.T.cell", "natural.killer.cell")
matching$B.plasma = c("B.cell","plasma.cell")
matching$mast = c("mast.cell")
matching$DC = c("dendritic.cell.type.1","dendritic.cell.type.2","activated.dendritic.cell","plasmacytoid.dendritic.cell")
matching$AEC = c("alveolar.epithelial.cell.type.1","alveolar.epithelial.cell.type.2")
#matching$other = c("blood.vessel.cell","ciliated.cell","muscle.cell","lymph.vessel.cell")
matching$muscle = c("muscle.cell")
matching$vessel = c("blood.vessel.cell", "lymph.vessel.cell")
matching$ciliated = c("ciliated.cell")
matching$fibroblast = c("fibroblast")


# actual decon
res <- spatialdecon(norm = as.matrix(assay(spe_keepNeg_ruv,"logcounts")),
                    bg = bg,
                    X = lung,
                    cellmerges = matching,
                    align_genes =TRUE)