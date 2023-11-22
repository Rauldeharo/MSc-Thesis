library(ggplot2)


#Load ChIP tables (fingeprint metrics with AUC value) ----
H3K4me1_chip <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/Fingerprint/metrics_H3K4me1_Fingerprint_ChIP.txt", sep = "\t", header = T)
H3K4me3_chip <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/Fingerprint/metrics_H3K4me3_Fingerprint_ChIP.txt", sep = "\t", header = T)
H3K27me3_chip <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/Fingerprint/metrics_H3K27me3_Fingerprint_ChIP.txt", sep = "\t", header = T)
H3K27ac_chip <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/Fingerprint/metrics_H3K27ac_Fingerprint_ChIP.txt", sep = "\t", header = T)
H2AK119ub_chip <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/Fingerprint/metrics_H2AK119ub_Fingerprint_ChIP.txt", sep = "\t", header = T)

#Creates dataframe with ChIP Histone mark and AUC value
chip_metrics <- data.frame()
chip_metrics <- rbind(H3K4me1_chip,H3K4me3_chip, H3K27me3_chip,H3K27ac_chip, H2AK119ub_chip )
chip_metrics <- aggregate(AUC ~ Sample, data = chip_metrics, FUN = mean) #to take the mean AUC between same inputs

histone <- substr(chip_metrics$Sample, 1, nchar(chip_metrics$Sample) - 2)

chip_metrics$Histone_Mark <- histone
chip_metrics <- data.frame(chip_metrics$Sample, chip_metrics$AUC, chip_metrics$Histone_Mark)
colnames(chip_metrics) <- c("Sample", "AUC","Histone_Mark")
chip_metrics$Omic <- "ChIP_seq"


chip_metrics$Histone_Mark <- factor(chip_metrics$Histone_Mark, levels = c("input","H3K4me1","H3K4me3","H3K27me3", "H3K27ac","H2AK119ub" ))

ggplot(data = chip_metrics,aes(x= Histone_Mark, y= AUC)) + geom_boxplot() + geom_jitter(aes(color = Histone_Mark),
                                                                                   width = 0.2,height = 0,
                                                                                   alpha = 0.7,
                                                                                   size = 3) +xlab("Histone Mark") +ylab("AUC") + theme(legend.position = "bottom")


#Load CnR tables (fingeprint metrics with AUC value) ----
H3K4me1_CnR <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/Fingerprint/metrics_H3K4me1_Fingerprint_CnR.txt", sep = "\t", header = T)
H3K4me3_CnR <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/Fingerprint/metrics_H3K4me3_Fingerprint_CnR.txt", sep = "\t", header = T)
H3K27me3_CnR <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/Fingerprint/metrics_H3K27me3_Fingerprint_CnR.txt", sep = "\t", header = T)
H3K27ac_CnR <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/Fingerprint/metrics_H3K27ac_Fingerprint_CnR.txt", sep = "\t", header = T)
H2AK119ub_CnR <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/Fingerprint/metrics_H2AK119ub_Fingerprint_CnR.txt", sep = "\t", header = T)

#Creates dataframe with CUT&RUN Histone mark and AUC value

CnR_metrics <- data.frame()
CnR_metrics <- rbind(H3K4me1_CnR,H3K4me3_CnR, H3K27me3_CnR,H3K27ac_CnR, H2AK119ub_CnR )
CnR_metrics <- aggregate(AUC ~ Sample, data = CnR_metrics, FUN = mean) #to take the mean between same inputs

histone <- substr(CnR_metrics$Sample, 1, nchar(CnR_metrics$Sample) - 2)

CnR_metrics$Histone_Mark <- histone
CnR_metrics <- data.frame(CnR_metrics$Sample, CnR_metrics$AUC, CnR_metrics$Histone_Mark)
colnames(CnR_metrics) <- c("Sample", "AUC","Histone_Mark")
CnR_metrics$Omic <- "CUT&RUN"


CnR_metrics$Histone_Mark <- factor(CnR_metrics$Histone_Mark, levels = c("input","H3K4me1","H3K4me3","H3K27me3", "H3K27ac","H2AK119ub" ))

ggplot(data = CnR_metrics,aes(x= Histone_Mark, y= AUC)) + geom_boxplot() + geom_jitter(aes(color = Histone_Mark),
                                                                                   width = 0.2,height = 0,
                                                                                   alpha = 0.7,
                                                                                   size = 3) +xlab("Histone Mark") +ylab("AUC") + theme(legend.position = "bottom")

#Global metrics (ChIP_seq and CUT&RUN) ----
global_metrics <- rbind(chip_metrics,CnR_metrics)

#Visualization of AUC values for both omics
ggplot(data = global_metrics,aes(x= Histone_Mark, y= AUC)) + geom_boxplot()+ facet_grid(. ~ Omic)+ geom_jitter(aes(color = Histone_Mark),
                    width = 0.2,height = 0,alpha = 0.7,size = 3) +xlab("Histone Mark") +ylab("AUC") + theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(color = "Histone mark") 





#Statistical comparison between input AUC values and Histone AUC values ----

##T-test for Chip ----
marks_chip <- unique(chip_metrics$Histone_Mark)
results_chip <- list()

for (mark in marks_chip) {
  input_aucs <- chip_metrics$AUC[chip_metrics$Histone_Mark == "input"]
  mark_aucs <- chip_metrics$AUC[chip_metrics$Histone_Mark == mark]
  
  t_test_result <- t.test(mark_aucs, input_aucs)
  results_chip[[mark]] <- t_test_result
  
}


##T-test for CUT&RUN ----
marks_CnR <- unique(CnR_metrics$Histone_Mark)
results_CnR <- list()

for (mark in marks_CnR) {
  input_aucs <- CnR_metrics$AUC[CnR_metrics$Histone_Mark == "input"]
  mark_aucs <- CnR_metrics$AUC[CnR_metrics$Histone_Mark == mark]
  
  t_test_result <- t.test(mark_aucs, input_aucs)
  results_CnR[[mark]] <- t_test_result
  
}

