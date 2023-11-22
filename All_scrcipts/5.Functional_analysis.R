library(GenomicRanges)
library(reshape2)
library(ggplot2)

#Function for better visualization
theme <- theme(text = element_text(size=11),
               strip.text.x = element_text(size = 11),
               axis.title=element_text(size=11),
               legend.text = element_text(size=11),
               axis.text =element_text(size=11),
               axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
               panel.background = element_rect(fill = "transparent"), # bg of the panel
               plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
               panel.grid.major = element_blank(), # get rid of major grid
               panel.grid.minor = element_blank(), # get rid of minor grid
               legend.background = element_rect(fill = "transparent"), # get rid of legend bg
               # legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
               strip.background =  element_rect(fill="transparent")
)

#Import promoters file and make GRanges object ----
promoter <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Files/promoter_regions_Homo_sapiens.GRCh38.104.bed", header = T)
GRpromoter <- makeGRangesFromDataFrame(promoter, keep.extra.columns = T)
strand(GRpromoter) <- "*"
seqlevelsStyle(GRpromoter) <- "UCSC"

#Import ChromHMM  output Bed files with states annotated ----

chip_10 <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/ChIP/nB_WT_10_dense_fin_chip.bed")
chip_10[is.na(chip_10)] <- "NA"
colnames(chip_10) <- c("seqnames", "start", "end", "state", "0", "strand")

CnR_10 <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/CnR/nB_WT_10_dense_fin.bed")
CnR_10[is.na(CnR_10)] <- "NA"
colnames(CnR_10) <- c("seqnames", "start", "end", "state", "0", "strand")

MACS_chip_10 <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/MACS2_ChIP/nB_WT_10_dense_fin.bed")
MACS_chip_10[is.na(MACS_chip_10)] <- "NA"
colnames(MACS_chip_10) <- c("seqnames", "start", "end", "state", "0", "strand")

MACS_CnR_10 <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/MACS2_CnR/nB_WT_10_dense_fin.bed")
MACS_CnR_10[is.na(MACS_CnR_10)] <- "NA"
colnames(MACS_CnR_10) <- c("seqnames", "start", "end", "state", "0", "strand")

SEACR_str_10 <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/SEACR_CnR_str/nB_WT_10_dense_fin.bed")
SEACR_str_10[is.na(SEACR_str_10)] <- "NA"
colnames(SEACR_str_10) <- c("seqnames", "start", "end", "state", "0", "strand")

SEACR_rel_10 <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/SEACR_CnR_rel/nB_WT_10_dense_fin.bed")
SEACR_rel_10[is.na(SEACR_rel_10)] <- "NA"
colnames(SEACR_rel_10) <- c("seqnames", "start", "end", "state", "0", "strand")



#Create list of States BED files
states_10 <- list(chip_10,CnR_10,MACS_chip_10,MACS_CnR_10,SEACR_str_10,SEACR_rel_10)
names(states_10) <- c("chip_10","CnR_10","MACS_chip_10","MACS_CnR_10","SEACR_str_10","SEACR_rel_10")

##Make Genomic Ranges from states BED files ----
GR_chip_10 <- makeGRangesFromDataFrame(chip_10)
GR_CnR_10 <- makeGRangesFromDataFrame(CnR_10)
GR_MACS_chip_10 <- makeGRangesFromDataFrame(MACS_chip_10)
GR_MACS_CnR_10 <- makeGRangesFromDataFrame(MACS_CnR_10)
GR_SEACR_str_10 <- makeGRangesFromDataFrame(SEACR_str_10)
GR_SEACR_rel_10 <- makeGRangesFromDataFrame(SEACR_rel_10)



#FIND INTERSECT REGIONS BETWEEN STATES REGIONS AND PROMOTERS REGIONS (to extract states-overlapping genes) ----
intersection <- list()
new_intersection <- list()

for (i in 1:length(states_10)) {
  GRanges_i <- makeGRangesFromDataFrame(states_10[[i]],keep.extra.columns = T)
  intersection[[i]] <- mergeByOverlaps(GRanges_i,GRpromoter,type = "within") #STATE REGIONS INSIDE THE PROMOTER REGIONS
  new_intersection[[i]] <- mergeByOverlaps(GRpromoter,GRanges_i,type = "within") #PROMOTER REGIONS INSIDE THE STATE REGIONS
  
}

names(intersection) <- paste0(names(states_10), "_inter")
names(new_intersection) <- paste0(names(states_10), "_inter")


for (i in 1:length(intersection)) {
  df1 <- as.data.frame(intersection[[i]])
  df2 <- as.data.frame(new_intersection[[i]])
  
  new_df1 <- df1[,c(1:6,16:25)]
  new_df1$condition <- "state_inside"
  new_df2 <- df2[,c(16:21,1:10)]
  new_df2$condition <- "promoter_inside"
  
  final_df <- rbind(new_df1,new_df2)
  write.table(x = final_df,file = paste0("~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/",names(intersection)[i],".tsv"),sep ="\t",col.names = T, row.names = F)
}
