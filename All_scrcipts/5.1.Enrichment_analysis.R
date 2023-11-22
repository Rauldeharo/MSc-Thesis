library(regioneReloaded)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)

#Load States-BED files ----
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

#Make GenomicRanges object: ----
GR_MACS_chip_10 <- makeGRangesFromDataFrame(MACS_chip_10,keep.extra.columns = T)
GR_MACS_chip_10 <- split(GR_MACS_chip_10,GR_MACS_chip_10$state)
names(GR_MACS_chip_10) <- paste(names(GR_MACS_chip_10),"MACS_ChIP", sep = "_")

GR_MACS_CnR_10 <- makeGRangesFromDataFrame(MACS_CnR_10,keep.extra.columns = T)
GR_MACS_CnR_10 <- split(GR_MACS_CnR_10,GR_MACS_CnR_10$state)
names(GR_MACS_CnR_10) <- paste(names(GR_MACS_CnR_10),"MACS_CnR", sep = "_")

GR_SEACR_str_10 <- makeGRangesFromDataFrame(SEACR_str_10,keep.extra.columns = T)
GR_SEACR_str_10 <- split(GR_SEACR_str_10,GR_SEACR_str_10$state)
names(GR_SEACR_str_10) <- paste(names(GR_SEACR_str_10),"SEACR_str", sep = "_")

GR_SEACR_rel_10 <- makeGRangesFromDataFrame(SEACR_rel_10,keep.extra.columns = T)
GR_SEACR_rel_10 <- split(GR_SEACR_rel_10,GR_SEACR_rel_10$state)
names(GR_SEACR_rel_10) <- paste(names(GR_SEACR_rel_10),"SEACR_rel", sep = "_")

#Load RNA-seq counts data: ----
counts <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Files/Counts_GeneName_norm_rlog_onlyBcells.tsv", header = T, sep = "\t")
nB_counts <- counts[,c("Geneid","nB_WT_1", "nB_WT_2", "nB_WT_3")] #keep only those counts in nB cells
nB_counts$mean <- rowMeans(nB_counts[,-1])

#Load Gene Promoters file annotated: ----
promoter <- read.table(file = "~/Documents/BIOINFORMATICS/TFM/Files/promoter_regions_Homo_sapiens.GRCh38.104.bed", header = T)
promoter$seqnames <- paste("chr",promoter$seqnames,sep = "")

##Take those genes present in Promoters file ----
nB_counts <- merge(nB_counts,promoter,by.x="Geneid", by.y="ID")

##Split genes based on their expression (quartiles), and make GRanges object ----
q_nB_WT <- quantile(nB_counts$mean)
q1_nB_WT <- makeGRangesFromDataFrame(nB_counts[nB_counts$mean <= q_nB_WT[2],c(1,6,7,8)])
q2_nB_WT <- makeGRangesFromDataFrame(nB_counts[(nB_counts$mean > q_nB_WT[2]) & (nB_counts$mean <= q_nB_WT[3]),c(1,6,7,8)])
q3_nB_WT <- makeGRangesFromDataFrame(nB_counts[(nB_counts$mean > q_nB_WT[3]) & (nB_counts$mean <= q_nB_WT[4]),c(1,6,7,8)])
q4_nB_WT <- makeGRangesFromDataFrame(nB_counts[(nB_counts$mean > q_nB_WT[4]) & (nB_counts$mean <= q_nB_WT[5]),c(1,6,7,8)])


#Create 2 sets for permutation test: ----

#First set: GRanges objects for the different States files
list_states <- list(GR_MACS_chip_10,GR_MACS_CnR_10, GR_SEACR_str_10, GR_SEACR_rel_10)
#Second set: Genes divided in 4 quartiles
list_genes <- list(q1=q1_nB_WT,q2=q2_nB_WT,q3=q3_nB_WT,q4=q4_nB_WT)



#Perform crosswise permutation test ----
enrich <-  crosswisePermTest(
  Alist = list_states,
  Blist = list_genes,
  min_sampling = 5000,
  ranFUN = "randomizeRegions", 
  evFUN = "numOverlaps", 
  ntimes = 1000,
  genome = "hg38",
  mc.cores = 20
)

##Generates the matrix summarizing the permutation test ----
enrich <- makeCrosswiseMatrix(enrich)

##Save enrichment plot ----
plot_enrich <- plotCrosswiseMatrix(enrich)

###Modifications for plot visualization ----
matrix_plot <- plot_enrich$data

matrix_plot$X <- factor(matrix_plot$X, levels = c("q1","q2", "q3", "q4"))

matrix_plot$condition <- gsub(pattern = ".*_",replacement = "",x = matrix_plot$Y)
matrix_plot$state <- gsub(pattern = "_.*",replacement = "",x = matrix_plot$Y)

matrix_plot$state<-factor(matrix_plot$state,levels = unique(matrix_plot$state)[c(3,1,2,9,12,13,6,4,10,5,8,11,7)])


ggplot(data = matrix_plot,aes(x= X, y = state, fill = value)) + 
  geom_tile() + 
  scale_fill_gradientn(colours=c("#ef6320","white","#032ead") ,
                       limits= c(-0.2,4),
                       values = scales::rescale(c(4,0,-0.2)), 
                       oob = scales::squish) +
  facet_wrap(~ condition, scales = "free_y")

