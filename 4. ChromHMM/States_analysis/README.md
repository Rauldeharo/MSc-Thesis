# Comparative analyisis between 5 to 10 states modes
With the aim to find out the most accurate number of states to be used in the LearnModel command, a comparative analysis between the results obtained from 5, 10 and 15 states was performed.

## Computing overlapped regionds between 5 and 10 states modes
Overlapped regions in the same state between 5 and 10 were computed.


```R
library(GenomicRanges)
library(ggplot2)

#Load States-BED files and preprocess for downstream analysis ----
bed5 <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/SEACR_CnR_str/nB_WT_5_dense_fin.bed")
bed10 <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/SEACR_CnR_str/nB_WT_10_dense_fin.bed")
bed15 <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/SEACR_CnR_str/nB_WT_15_dense_fin.bed")

#Replace NA (missing values) for "NA" state value
bed5$V4[is.na(bed5$V4)] <- "NA"
bed10$V4[is.na(bed10$V4)] <- "NA"
bed15$V4[is.na(bed15$V4)] <- "NA"

#Naming columns of BED files to further create GRanges
names(bed5) <- c("chr", "start", "end", "state", "0", "strand", "met_start", "met_end", "color")
names(bed10) <- c("chr", "start", "end", "state", "0", "strand", "met_start", "met_end", "color")
names(bed15) <- c("chr", "start", "end", "state","0", "strand", "met_start", "met_end", "color")

#Substract 1 base at the end of each region
bed5$end <- bed5$end -1
bed10$end <- bed10$end -1
bed15$end <- bed15$end -1



#Find overlapping regions in the same state between 5state and 10state: ----
overlap_table <- data.frame(state = character(),
                            coverage_inter = numeric(), 
                            coverage_GR5 = numeric(),
                            coverage_GR10 = numeric(),
                            percent_inter_to10 = numeric(),
                            percent_inter_to5 = numeric(),
                            percent_new_states_in10 = numeric(),
                            stringsAsFactors = FALSE)

for (state in unique(bed5$state)) {
  
  #Subset bed files for the current state
  bed5_state <- bed5[bed5$state == state,]
  bed10_state <- bed10[bed10$state == state,]

  #Create GRanges from bed files in each state
  GR5 <- makeGRangesFromDataFrame(bed5_state,keep.extra.columns = TRUE)
  GR10 <- makeGRangesFromDataFrame(bed10_state,keep.extra.columns = TRUE)
  
  #Compute coverage for each Granges
  bed_intersect <- intersect(GR5,GR10)

  coverage_inter <- sum(width(bed_intersect))

  coverage_GR5 <- sum(width(GR5))

  coverage_GR10 <- sum(width(GR10))
  
  percent_inter_to5 <- (coverage_inter/coverage_GR5)*100
  percent_inter_to10 <- (coverage_inter/coverage_GR10)*100
  
  percent_new_states_in10 <- (100-percent_inter_to5)  
  
  results <- c(state,coverage_GR5, coverage_GR10, coverage_inter ,percent_inter_to5 ,percent_inter_to10, percent_new_states_in10)
  
  #Add results in the overlap table (row by row)
  overlap_table <- rbind(overlap_table,results)
  colnames(overlap_table) <- c("state",
                               "coverage_GR5", "coverage_GR10", "coverage_inter","percent_inter_to5","percent_inter_to10",
                               "percent_new_states_in10")
  }


write.table(overlap_table, file = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/SEACR_CnR_str/overlap_SEACR_CnR_str_5_10.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



```

## New states

Compute which states appears as new in 10 states mode that not appears in 5 states mode.

```R
#Find new states in 10states that not appears in 5states: ----

#Create global Granges for 10 states
GR10_glob <- makeGRangesFromDataFrame(bed10,keep.extra.columns = TRUE)

df_final <- data.frame()

for (state in unique(bed5$state)) {
  
  #Subset bed5 for the current state and creates GRanges for this bed5_state
  bed5_state <- bed5[bed5$state == state,]
  GR5 <- makeGRangesFromDataFrame(bed5_state,keep.extra.columns = TRUE)
  
  #Find disjoint regions between GR5, GR10_glob
  dis1 <- disjoin(c(GR5,GR10_glob), with.revmap=T)
  dis <- dis1[lapply(dis1$revmap,length) > 1] #takes those that comes from both Granges
  bed_intersect <- mergeByOverlaps(dis, GR10_glob)
 
  
  df <- as.data.frame(table(unique(bed_intersect$state)))

  for (st in unique(bed_intersect$state)){
    
    cov_bed_intersect <- subset(bed_intersect, state == st)
    coverage_inter <- sum(width(cov_bed_intersect$dis))
    df$coverage[df$Var1 == st] <- coverage_inter
    
  }
  
  df$state5 <- state
  df_final <- rbind(df_final,df)

}
colnames(df_final) <- c("state10", "Num.regions_in10", "Coverage_Regions_in10","state5")

write.table(df_final, file = "~/Documents/BIOINFORMATICS/TFM/Results/ChromHMM/SEACR_CnR_str/new_states_SEACR_CnR_str_5_10.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

ggplot(df_final,aes(x=state5,y=Coverage_Regions_in10,fill=state10)) + 
  geom_bar(stat="identity",position="fill") + ggtitle("New states in 10 over 5") +
  ggpubr::theme_classic2()
