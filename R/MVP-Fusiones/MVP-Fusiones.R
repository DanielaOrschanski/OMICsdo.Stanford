#Librerias:
library(shiny)
library(shinydashboard)
library(DT)
library(readxl)
library(dplyr)
#install.packages("shinyFiles")
library(shinyFiles)
#install.packages("shinycssloaders")
library(shinycssloaders)
library(plotly)

#Todos los reportes de fusiones:
ChinaSRA <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Todos-FusionReports_ChinaSRA.xlsx")
colnames(ChinaSRA)[1] <- "ID"
ChinaSRA$Cohort <- "Breast_China"
ChinaSRA$Cancer <- "Breast"

Gastrico_5p <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Todos-FusionReports_Gastrico-5p.xlsx")
colnames(Gastrico_5p)[1] <- "ID"
Gastrico_5p$Cohort <- "Gastric_5p"
Gastrico_5p$Cancer <- "Gastric"

Gastrico_6p <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Todos-FusionReports_Gastrico-6p.xlsx")
colnames(Gastrico_6p)[1] <- "ID"
Gastrico_6p$Cohort <- "Gastric_6p"
Gastrico_6p$Cancer <- "Gastric"

Gastrico_610 <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Todos-FusionReports_Gastrico-610.xlsx")
colnames(Gastrico_610)[1] <- "ID"
Gastrico_610$Cohort <- "Gastric_610"
Gastrico_610$Cancer <- "Gastric"
Metadata_Gastrico_610 <- read_csv("/media/16TBDisk/Daniela/Fusiones/Gastrico-610/Metadata-Gastrico-610.csv")
lineas_celulares <- Metadata_Gastrico_610$Run[which(Metadata_Gastrico_610$source_name == "Gastric cancer cell line")]
Gastrico_610 <- Gastrico_610[-which(Gastrico_610$ID %in% lineas_celulares),]

library(OMICsdo)
library(openxlsx)
#QC_seq_align(patients_dir <-  "/media/16TBDisk/Daniela/Fusiones/Gastrico-610")


Leucemia_1_ <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Todos-FusionReports_Leucemia.xlsx")
colnames(Leucemia_1_)[1] <- "ID"
Leucemia_1_$Cohort <- "Leukemia-Pat"
Leucemia_1_$Cancer <- "Leukemia"

Tiroides_Agresivo_446 <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Todos-FusionReports_Tiroides-446.xlsx")
colnames(Tiroides_Agresivo_446)[1] <- "ID"
Tiroides_Agresivo_446$Cohort <- "AgressiveThyroid_446"
Tiroides_Agresivo_446$Cancer <- "Thyroid"


Tiroides_PRJNA695543 <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Todos-FusionReports_Tiroides-PRJNA695543.xlsx")
colnames(Tiroides_PRJNA695543)[1] <- "ID"
Tiroides_PRJNA695543$Cohort <- "Thyroid_543"
Tiroides_PRJNA695543$Cancer <- "Thyroid"
Tiroides_PRJNA695543 <- Tiroides_PRJNA695543[, -which(colnames(Tiroides_PRJNA695543) %in% c("MTT", "Grupo"))]

Melanoma <-read_csv("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Todos-FusionReports_Melanoma.csv")
colnames(Melanoma)[1] <- "ID"
Melanoma$Cohort <- "Melanoma"
Melanoma$Cancer <- "Melanoma"
Melanoma <- Melanoma[, -which(colnames(Melanoma) %in% c("Fus", "group"))]


FLENI <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Todos-FusionReports_FLENI.xlsx")
colnames(FLENI)[1] <- "ID"
FLENI$Cohort <- "FLENI"
FLENI$Cancer <- "CNS"

Gastrico80 <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Todos-FusionReports_Gastric-80p.xlsx")
colnames(Gastrico80)[1] <- "ID"
Gastrico80$Cohort <- "Gastric-80p"
Gastrico80$Cancer <- "Gastric"

Tir591 <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Todos-FusionReports_Tiroides-PRJEB11591.xlsx")
colnames(Tir591)[1] <- "ID"
Tir591$Cohort <- "Thyroid-591"
Tir591$Cancer <- "Thyroid"

ReportFusiones_TodasCohortes <- rbind(FLENI, Tiroides_Agresivo_446, Tiroides_PRJNA695543,
                                      Leucemia_1_, Gastrico_5p, Gastrico_6p, Gastrico_610,
                                      ChinaSRA, Melanoma, Gastrico80, Tir591)


library(openxlsx)
write.xlsx(ReportFusiones_TodasCohortes, file = "/media/16TBDisk/Daniela/Fusiones/Todos_Reports/ReportFusions_TodasCohortes.xlsx")
write.xlsx(ReportFusiones_TodasCohortes, file = "/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/Act_ReportFusions_TodasCohortes.xlsx")

ReportFusiones_TodasCohortes <- read.xlsx("/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/Act_ReportFusions_TodasCohortes.xlsx")
length(unique(ReportFusiones_TodasCohortes$Cohort))


################################################################################

GeneFusions_Gastrico_610 <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Genes/GeneFusions_Gastrico-610.xlsx")
unique(ReportFusiones_TodasCohortes$Cohort)
GeneFusions_Gastrico_610$Porcentaje_Apariciones <- (GeneFusions_Gastrico_610$Total_Apariciones / nrow(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "Gastric_610"),]))*100
GeneFusions_Gastrico_610$Porcentaje_Muestras_Distintas <- (GeneFusions_Gastrico_610$Muestras_Distintas / length(unique(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "Gastric_610"),]$ID)))*100
GeneFusions_Gastrico_610$Cohort <- "Gastric_610"

GeneFusions_China <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Genes/GeneFusions_China.xlsx")
unique(ReportFusiones_TodasCohortes$Cohort)
GeneFusions_China$Porcentaje_Apariciones <- (GeneFusions_China$Total_Apariciones / nrow(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "Breast_China"),]))*100
GeneFusions_China$Porcentaje_Muestras_Distintas <- (GeneFusions_China$Muestras_Distintas / length(unique(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "Breast_China"),]$ID)))*100
colnames(GeneFusions_China)
GeneFusions_China <-GeneFusions_China[, c(1,2,3,11,12)]
GeneFusions_China$Cohort <- "Breast_China"

GeneFusions_Tir446 <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Genes/GeneFusions_Tiroides-446.xlsx")
unique(ReportFusiones_TodasCohortes$Cohort)
GeneFusions_Tir446$Porcentaje_Apariciones <- (GeneFusions_Tir446$Total_Apariciones / nrow(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "AgresiveThyroid_446"),]))*100
GeneFusions_Tir446$Porcentaje_Muestras_Distintas <- (GeneFusions_Tir446$Muestras_Distintas / length(unique(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "AgresiveThyroid_446"),]$ID)))*100
GeneFusions_Tir446$Cohort <- "AgressiveThyroid_446"


GeneFusions_G80 <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Genes/GeneFusions_Gastric-80p.xlsx")
unique(ReportFusiones_TodasCohortes$Cohort)
colnames(GeneFusions_G80)[2] <- "Total_Apariciones"
colnames(GeneFusions_G80)[3] <- "Muestras_Distintas"
GeneFusions_G80$Porcentaje_Apariciones <- (GeneFusions_G80$Total_Apariciones / nrow(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "Gastric-80p"),]))*100
GeneFusions_G80$Porcentaje_Muestras_Distintas <- (GeneFusions_G80$Muestras_Distintas / length(unique(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "Gastric-80p"),]$ID)))*100
GeneFusions_G80$Cohort <- "Gastric-80p"


GeneFusions_FLENI <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Genes/GeneFusions_FLENI.xlsx")
unique(ReportFusiones_TodasCohortes$Cohort)
GeneFusions_FLENI$Porcentaje_Apariciones <- (GeneFusions_FLENI$Total_Apariciones / nrow(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "FLENI"),]))*100
GeneFusions_FLENI$Porcentaje_Muestras_Distintas <- (GeneFusions_FLENI$Muestras_Distintas / length(unique(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "FLENI"),]$ID)))*100
GeneFusions_FLENI$Cohort <- "FLENI"

GeneFusions_Tir591 <- read_excel("/media/16TBDisk/Daniela/Fusiones/Todos_Reports/Genes/GeneFusions_Tiroides-PRJEB11591.xlsx")
colnames(GeneFusions_Tir591)[2] <- "Total_Apariciones"
colnames(GeneFusions_Tir591)[3] <- "Muestras_Distintas"
GeneFusions_Tir591$Porcentaje_Apariciones <- (GeneFusions_Tir591$Total_Apariciones / nrow(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "Thyroid-591"),]))*100
GeneFusions_Tir591$Porcentaje_Muestras_Distintas <- (GeneFusions_Tir591$Muestras_Distintas / length(unique(ReportFusiones_TodasCohortes[which(ReportFusiones_TodasCohortes$Cohort == "Thyroid-591"),]$ID)))*100
GeneFusions_Tir591$Cohort <- "Thyroid-591"

Genes_TodasCohortes <- rbind(GeneFusions_G80, GeneFusions_Tir446, GeneFusions_China,
                             GeneFusions_Gastrico_610, GeneFusions_FLENI, GeneFusions_Tir591)
colnames(Genes_TodasCohortes)[2] <- "Fused_Frequency"
colnames(Genes_TodasCohortes)[3] <- "Fused_Frequency_Per_ID"
colnames(Genes_TodasCohortes)[4] <- "Fused_Frequency_%"
colnames(Genes_TodasCohortes)[5] <- "Fused_Frequency_Per_ID_%"

Genes_TodasCohortes$`Fused_Frequency_%` <- round(Genes_TodasCohortes$`Fused_Frequency_%`,2)
Genes_TodasCohortes$`Fused_Frequency_Per_ID_%` <- round(Genes_TodasCohortes$`Fused_Frequency_Per_ID_%`,2)

library(openxlsx)
write.xlsx(Genes_TodasCohortes, file = "/media/16TBDisk/Daniela/Fusiones/Todos_Reports/ReportGenes_TodasCohortes.xlsx")
write.xlsx(Genes_TodasCohortes, file = "/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/ReportGenes_TodasCohortes.xlsx")

Genes_TodasCohortes <-  read_excel("/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/ReportGenes_TodasCohortes.xlsx")

############################################################################


