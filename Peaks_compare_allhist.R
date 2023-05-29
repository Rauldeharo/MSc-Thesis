library(GenomicRanges)
library(ggplot2)
library(reshape2)

## Functions----
compute_overlap <- function(first,second) #returns a list (row) with % of interaction (Xi,Y1) (Xi,y2) (Xi,y3)...
{
  sum_x <- sum(width(first))
  int <- GenomicRanges::intersect(first,second)
  overlap <- sum(width(int))*100/sum_x
  return(overlap)
}

compute_jaccard <- function(first,second)
{
  sum_x <- sum(width(first))
  sum_y <- sum(width(second))
  sum_int <- sum(width(intersect(first,second)))
  jac_idx <- (sum_int /((sum_x + sum_y)-sum_int))*100
  return(jac_idx)
}

compute_ov.coefficient <- function(first,second)
{
  sum_x <- sum(width(first))
  sum_y <- sum(width(second))
  sum_int <- sum(width(intersect(first,second)))
  if (sum_x < sum_y){
    ov.coeff <- (sum_int/sum_x)*100
    return(ov.coeff)
  } else {
    ov.coeff <- (sum_int/sum_y)*100
    return(ov.coeff)
  }
}



setwd("~/Documents/BIOINFORMATICS/TFM/Files/Histone_marks/")

##Create a vector of files----
files<- list.files("~/Documents/BIOINFORMATICS/TFM/Files/Histone_marks/", full.names = TRUE)

plots.similarity <- list() #Creates a list to store similarity plots
plots.jaccard <- list() #Creates a list to store jaccard plots
plots.coefficient <- list() #Creates a list to store ov.coefficient plots

overlaps_list <- list() #Creates a list to store overlaps_Hist_mark data.frames
jaccard_hist_list <- list() #Creates a list to store jaccard index data.frames
coeff_hist_list <- list() #Creates a list to store ov.coefficient index data.frames

##Creates a vector with the name of each histone mark----
for (mark in c("H3K4me3","H3K4me1","H3K27me3","H3K27ac","H2AK119ub")){
  sub_files <- grep(mark,files,value = T)
  sub_files<- sub_files[order(reverse(sub_files),decreasing = TRUE)]

  Hist_mark <- list() #Empty List for each Histone mark
  
  summary_peaks <- data.frame(col1 = character(),col2 = numeric(), col3 = numeric(),stringsAsFactors = FALSE)
  measures <- data.frame(col1 = character(),col2 = numeric(), col3 = numeric(),col4 = numeric(),
                            col5 = numeric(), col6 = numeric(), col7 = numeric(),stringsAsFactors = FALSE)
  
  #Iterate over histone mark: comparison between other peak calling or replicates---- 
  for (file in sub_files){
    x_df <- read.table(file,header = FALSE, sep="\t",stringsAsFactors = FALSE)
    x_gr <- makeGRangesFromDataFrame(x_df, keep.extra.columns = TRUE, seqnames.field = "V1", start.field = "V2", end.field = "V3")
    Hist_mark[[basename(file)]] <- x_gr #to store the GRanges object in Hist_mark GRangesList
    summary_peaks <- rbind(summary_peaks,c(basename(file),length(x_gr@ranges),sum(x_gr@ranges@width)))
    
    stats <- summary(x_gr@ranges@width) #width peak average
    measures <- rbind(measures, c(basename(file), stats[["Min."]], stats[["Max."]], stats[["1st Qu."]],stats[["Median"]], stats[["3rd Qu."]], stats[["Mean"]]))
  }
  colnames(summary_peaks) <- c("Sample","Peaks","Coverage")
  colnames(measures) <- c("Sample", "Min.","Max.","1st Qu.","Median", "3rd Qu.", "Mean")
  
  
  statistics <- merge(summary_peaks,measures, all=T)
  
  statistics$Sample <- gsub("nB_WT_","",statistics$Sample)
  statistics$Sample <- gsub("\\.bed|\\.narrow.*|\\.broad.*","",statistics$Sample)
  
  write.table(statistics,file = paste0("~/Documents/BIOINFORMATICS/TFM/Results/Metrics/statistics_",mark,".tsv"),
              row.names = F, col.names = T, quote = F, sep = "\t")
  
  #Apply function
  l <- lapply(Hist_mark, function(x)
  {
    ov <- lapply(Hist_mark, function(y) compute_overlap(x,y))
    as.data.frame(ov)
  
  })
  overlaps_Hist_mark <- as.data.frame(data.table::rbindlist(l)) #merge all data.frames (row) in one single data.frame
  rownames(overlaps_Hist_mark) <- colnames(overlaps_Hist_mark)
 
  m <- lapply(Hist_mark, function(x)
  {
    jc <- lapply(Hist_mark, function(y) compute_jaccard(x,y))
    as.data.frame(jc)
  })
  jaccard_hist <- as.data.frame(data.table::rbindlist(m))
  rownames(jaccard_hist) <- colnames(jaccard_hist)
  
  
  n <- lapply(Hist_mark,function(x)
    {
    coeff <- lapply(Hist_mark,function(y) compute_ov.coefficient(x,y))
    as.data.frame(coeff)
  })
  coeff_hist <- as.data.frame(data.table::rbindlist(n))
  rownames(coeff_hist) <- colnames(coeff_hist)

  ##Modifications for visualization----
  #Overlap:
  colnames(overlaps_Hist_mark)<- gsub("nB_WT_","",colnames(overlaps_Hist_mark))
  colnames(overlaps_Hist_mark)<- gsub("\\.bed|\\.narrow.*|\\.broad.*","",colnames(overlaps_Hist_mark))
  
  rownames(overlaps_Hist_mark)<- gsub("nB_WT_","",rownames(overlaps_Hist_mark))
  rownames(overlaps_Hist_mark)<- gsub("\\.bed|\\.narrow.*|\\.broad.*","",rownames(overlaps_Hist_mark))
  
  #Jaccard
  colnames(jaccard_hist)<- gsub("nB_WT_","",colnames(jaccard_hist))
  colnames(jaccard_hist)<- gsub("\\.bed|\\.narrow.*|\\.broad.*","",colnames(jaccard_hist))
  
  rownames(jaccard_hist)<- gsub("nB_WT_","",rownames(jaccard_hist))
  rownames(jaccard_hist)<- gsub("\\.bed|\\.narrow.*|\\.broad.*","",rownames(jaccard_hist))
  
  #Ov.coefficient
  colnames(coeff_hist)<- gsub("nB_WT_","",colnames(coeff_hist))
  colnames(coeff_hist)<- gsub("\\.bed|\\.narrow.*|\\.broad.*","",colnames(coeff_hist))
  
  rownames(coeff_hist)<- gsub("nB_WT_","",rownames(coeff_hist))
  rownames(coeff_hist)<- gsub("\\.bed|\\.narrow.*|\\.broad.*","",rownames(coeff_hist))

  
  #Transform wide dataframe to long dataframe (to visualize it by ggplot):
  overlaps_Hist_mark$row <- names(overlaps_Hist_mark) #create new variable row (with row names) --> needed for melt function
  overlaps_list[[mark]] <- overlaps_Hist_mark #list of data.frames (overlaps_Hist_mark)

  jaccard_hist$row <- names(jaccard_hist)
  jaccard_hist_list[[mark]] <- jaccard_hist #list of data.frames (jaccard)
  
  coeff_hist$row <- names(coeff_hist)
  coeff_hist_list[[mark]] <- coeff_hist #list of data.frames (jaccard)
  
  #Creates tables (overlap_Hist_mark, jaccard index, ov.coefficient) in my directory----
  write.table(overlaps_Hist_mark,file = paste0("~/Documents/BIOINFORMATICS/TFM/Results/Metrics/overlap_",mark,".tsv"),
              row.names = F, col.names = T, quote = F, sep = "\t")
  
  write.table(jaccard_hist,file = paste0("~/Documents/BIOINFORMATICS/TFM/Results/Metrics/jaccard_",mark,".tsv"),
              row.names = F, col.names = T, quote = F, sep = "\t")
  
  write.table(coeff_hist,file = paste0("~/Documents/BIOINFORMATICS/TFM/Results/Metrics/ov_coefficient_",mark,".tsv"),
              row.names = F, col.names = T, quote = F, sep = "\t")
  
  long_ov <- melt(overlaps_Hist_mark) #use the only variable that is a factor (row variable)
  names(long_ov) <- c("P1","P2","Sim") #rename the long_ov dataframe recently created
  
  long_jac <- melt(jaccard_hist)
  names(long_jac) <- c("P1","P2","Jac_idx")
  
  long_ov.coef <- melt(coeff_hist)
  names(long_ov.coef) <- c("P1","P2","Ov.coeff")
  
  #Creates inner order (factor = levels):
  long_ov$P1 <- factor(long_ov$P1,levels = unique(long_ov$P1)[c(3,2,1,6,5,4,12,11,10,9,8,7)])
  long_ov$P2 <- factor(long_ov$P2,levels = unique(long_ov$P2)[c(3,2,1,6,5,4,12,11,10,9,8,7)])
  long_ov$Sim <- as.character((long_ov$Sim))
  
  long_jac$P1 <-  factor(long_jac$P1,levels = unique(long_jac$P1)[c(3,2,1,6,5,4,12,11,10,9,8,7)])
  long_jac$P2 <- factor(long_jac$P2,levels = unique(long_jac$P2)[c(3,2,1,6,5,4,12,11,10,9,8,7)])
  long_jac$Jac_idx <- as.character(long_jac$Jac_idx)

  long_ov.coef$P1 <-  factor(long_ov.coef$P1,levels = unique(long_ov.coef$P1)[c(3,2,1,6,5,4,12,11,10,9,8,7)])
  long_ov.coef$P2 <- factor(long_ov.coef$P2,levels = unique(long_ov.coef$P2)[c(3,2,1,6,5,4,12,11,10,9,8,7)])
  long_ov.coef$Ov.coeff <- as.character(long_ov.coef$Ov.coeff, digits = 2)
  
  
  #Visualization: ggplot (introducing  data.frames with inner factor order)----
  similarity <- ggplot(data = long_ov, aes(x=P2,y=P1)) + geom_tile(aes(fill = as.numeric(Sim))) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1, size = 3 ),axis.text.y = (element_text(size = 3)), 
          legend.title = element_text(size = 4),
          legend.text = element_text(size = 4),axis.title.x = element_blank(),
          axis.title.y = element_blank(),plot.title = element_text(size = 6)) + geom_text(aes(label = round(as.numeric(Sim),digits = 1)),position = position_dodge(width=0.3),  size=1.5) +
    scale_fill_gradientn(colours = c("#314d1d","#d9fb83","white"),limits= c(0,100),values = scales::rescale(c(100,50,0)),oob = scales::squish ,name= "Overlap %") +labs(title =paste0(mark," Overlap"))
  
  plots.similarity[[mark]] <- similarity

  jaccard_idx <- ggplot(data = long_jac, aes(x=P2,y=P1)) + geom_tile(aes(fill = as.numeric(Jac_idx))) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1, size = 3 ),axis.text.y = (element_text(size = 3)), 
          legend.title = element_text(size = 4),
          legend.text = element_text(size = 4),axis.title.x = element_blank(),
          axis.title.y = element_blank(),plot.title = element_text(size = 6)) + geom_text(aes(label = round(as.numeric(Jac_idx),digits = 1)),position = position_dodge(width=0.3),  size=1.5)+
    scale_fill_gradientn(colours = c("#0D47A1","#BBDEFB","white"),limits= c(0,100),values = scales::rescale(c(100,50,0)),oob = scales::squish ,name= "Jaccard index %")+labs(title = paste0(mark," Jaccard"))
  
  plots.jaccard[[mark]] <- jaccard_idx 
  
  
  ov.coefficient <- ggplot(data = long_ov.coef, aes(x=P2,y=P1)) + geom_tile(aes(fill = as.numeric(Ov.coeff))) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1, size = 3 ),axis.text.y = (element_text(size = 3)), 
          legend.title = element_text(size = 4),
          legend.text = element_text(size = 4),axis.title.x = element_blank(),
          axis.title.y = element_blank(),plot.title = element_text(size = 6)) + geom_text(aes(label = round(as.numeric(Ov.coeff),digits = 1)),position = position_dodge(width=0.3),  size=1.5)+
    scale_fill_gradientn(colours=c("#a88102","#ffe9a1","white") ,limits= c(0,100),values = scales::rescale(c(100,50,0)), oob = scales::squish, name = "Overlap coefficient %" )+labs(title = paste0(mark," Overlap coefficient"))
  
  plots.coefficient[[mark]] <- ov.coefficient
  
}

#Plotting in different pdf files:
pdf("~/Documents/BIOINFORMATICS/TFM/Results/Metrics/overlap_plots.pdf",height = 3,width = 4)
for (i in plots.similarity){
  print(i)
}
dev.off()

pdf("~/Documents/BIOINFORMATICS/TFM/Results/Metrics/jaccard_plots.pdf",height = 3,width = 4)
for (j in plots.jaccard){
  print(j)
}
dev.off()

pdf("~/Documents/BIOINFORMATICS/TFM/Results/Metrics/ov_coefficient.pdf",height = 3,width = 4)
for (k in plots.coefficient){
  print(k)
}
dev.off()


