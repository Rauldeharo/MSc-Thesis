#!/bin/bash

#Fingerprint for H3K4me1 performed by CUT&RUN


#Histone mark to analyze
Histone_mark="H3K4me1"
plot_title=$Histone_mark"_Fingerprint_CnR"

#Paths to BAM files (samples and input/IgG controls)
path_samples="/mnt/beegfs/rdeharo/cutnrun_link_samples/"
path_controls="/mnt/beegfs/rdeharo/cutnrun_link_controls/"

#change file names for each histone mark
sampleA="nB_WT_H3K4me1_2_CUTnRUN.filt.nodup.bam"
sampleB="nB_WT_H3K4me1_3_CUTnRUN.filt.nodup.bam"
sampleC="nB_WT_H3K4me1_8_CUTnRUN.filt.nodup.bam"

controlA="nB_WT_IgG_2_CUTnRUN.filt.nodup.bam"
controlB="nB_WT_IgG_3_CUTnRUN.filt.nodup.bam"
controlC="nB_WT_IgG_8_CUTnRUN.filt.nodup.bam"

#change replicate number for each histone mark
label_sampA=$Histone_mark"_2"
label_sampB=$Histone_mark"_3"
label_sampC=$Histone_mark"_8"

#change input number for each histone mark
label_inputA="input_2"
label_inputB="input_3"
label_inputC="input_8"

plotFingerprint -b $path_samples$sampleA $path_samples$sampleB $path_samples$sampleC $path_controls$controlA $path_controls$controlB $path_controls$controlC -T $plot_title --plotFile /mnt/beegfs/rdeharo/Fingerprint/$plot_title.pdf --outQualityMetrics --outRawCounts /mnt/beegfs/rdeharo/Fingerprint/$plot_title.tab --labels $label_sampA $label_sampB $label_sampC $label_inputA $label_inputB $label_inputC
