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



