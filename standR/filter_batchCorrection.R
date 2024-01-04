### remove ROIs with low % alignement OR sequencing saturation

qc <- colData(speGQC)$AlignedReads/colData(speGQC)$RawReads >=0.9 & colData(speGQC)$SequencingSaturation >=90
table(qc)
speFinal <- speGQC[,qc]

### batch correction with RUV
speRUV <- findNCGs(speFinal, batch_name="patid", top_n=200)
for(i in seq(4)){
  spe_ruv <- geomxBatchCorrection(speRUV, factors = "type", 
                                  NCGs = metadata(speRUV)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = patid, title = paste0("k = ", i)))
  
}


spe_ruv <- geomxBatchCorrection(speRUV, factors = "anno_type", 
                                NCGs = metadata(speRUV)$NCGs, k = 4)


saveRDS(spe_ruv, "spe_ruv.rds")
