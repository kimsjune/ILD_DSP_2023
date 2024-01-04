# importing RAW tables downloaded from GeoMx 
# R packages cannot generate the same tables as GeoMx website can.

library(standR)
setwd("C:/Users/jk/Documents/Lab2/nanostr/dge_final")

### reading tables
# copy from Excel into Notepad++ to avoid an EOL error
segProp <- read.table("./../dsp_download/segProp_remove_GLI_with_outcome.txt",
                      header=F,
                      sep="\t")

bioprob <- read.table("./../dsp_download/bioprobe_remove_GLI.txt",
                      header=F,
                      sep="\t")

countFile <- as.data.frame(bioprob[,c(3,13:191)])
sampleAnnoFile <- as.data.frame(segProp)
featureAnnoFile <- as.data.frame( bioprob[,c(1:12)])


# for some reason, thes tables must be written out as a file then imported directly into readGeoMx
write.table(countFile,file="export_countFile.txt",sep="\t",
            row.names=F, col.names=F)
write.table(sampleAnnoFile,file="export_sampleAnnoFile.txt",sep="\t",row.names=F, col.names=F)
write.table(featureAnnoFile,file="export_featureAnnoFile.txt",sep="\t",row.names=F, col.names=F)

spe <- readGeoMx("./export_countFile.txt",
                 "./export_sampleAnnoFile.txt",
                 "./export_featureAnnoFile.txt")
