# Functional analysis integrating ChromHMM states with RNA-seq
In order to find a biological meaning in the states annotated as active or repressed, merge of ChromHMM states previously obtained with RNA-seq data was performed. Therefore, a states-overlapping gene promoters (defined as -1000 +250bp from the TSS) were used to determine their correlation between state classifcation and RNA-seq expression.

## Grouping states-overlapping gene promoters
First, genes (annotated based on promoter regions) that overlap state regions were selected and classified into active-repressed genes (based on the overlapped state)
```R
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
```

## Merging with RNA-seq data
Then, genes previously selected were correlated with their expression in the RNA-seq.
```R
#MERGE WITH RNA-seq ----
library(tidyverse)

##import RNA-seq counts table ----
counts <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Files/Counts_GeneName_norm_rlog_onlyBcells.tsv", header = T, sep = "\t")
nB_counts <- counts[,c("Geneid","nB_WT_1", "nB_WT_2", "nB_WT_3")] #keep only those counts in nB cells

##define active and repressive states ----
repressive_states <- c("PCdom", "PCdomA", "PCdomBiv", "ReprPC", "ReprPC2", "NA")
active_states <- c("TssA", "TssFlnk", "TssFlnkU","TssFlnkD")

##import recently created files State regions - Promoter regions ----
inter_files <- list.files(path = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM",pattern = "\\.tsv$", full.names = TRUE)

all_conditions <- data.frame()

for (inter in inter_files) {
  genes_10_states <- read.table(file = inter, header = T, sep = "\t")
  genes_10_states <- as_tibble(genes_10_states)
  genes_10_states<- genes_10_states %>% group_by(GRpromoter.ID) %>% summarise(stateSum=paste(GRanges_i.state,collapse = ","))
  
  genes_10_states <- subset(genes_10_states,genes_10_states$stateSum != "NA")
  
  active_genes_10_states <- subset(genes_10_states, unlist(lapply(genes_10_states$stateSum, function(x) any(str_split(x,",")[[1]] %in% active_states))))
  active_genes_10_states <- subset(nB_counts, Geneid %in% active_genes_10_states$GRpromoter.ID)

  if (nrow(active_genes_10_states) != 0) #avoiding errors in files that not contain active states 
  {
    active_genes_10_states$State <- "Active"
  }
  
  repressive_genes_10_states <- subset(genes_10_states, unlist(lapply(genes_10_states$stateSum, function(x) all(str_split(x,",")[[1]] %in% repressive_states))))
  repressive_genes_10_states <- subset(nB_counts, Geneid %in% repressive_genes_10_states$GRpromoter.ID)
  repressive_genes_10_states$State <- "Repressive"
  
  final_genes_10_states <- rbind(active_genes_10_states,repressive_genes_10_states)
  
  final_genes_10_states <- melt(final_genes_10_states)
  final_genes_10_states$value1 <- final_genes_10_states$value +1
  
  df_name<- basename(inter)
  assign(df_name, final_genes_10_states)
  
  
  cond <- gsub("_10_inter.tsv","",df_name)
  final_genes_10_states$condition <- cond
  
  all_conditions <- rbind(all_conditions,final_genes_10_states)
  
  
}

#Visualization of correlation between states and gene expression ----
ggplot(data = all_conditions, aes(x = variable,y = value, fill = State))+ geom_boxplot() +ylab("log10 mRNA expression")+ facet_wrap(~condition)+theme+ theme(panel.border=element_rect(color="black",fill=NA))
```
