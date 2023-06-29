# File modification with the purpose to use SEACR peak-calling:


## Transform BAM to a Bedgraph file (input to SEACR peak calling)
From BAM filtered file (same pipeline as Chip-seq) we get a a BAM No duplicates file.
We have to modify it to use bedraph file in SEACR peak-caller

``` !/bin/bash

input_file="/mnt/beegfs/rdeharo/prueba/nB_WT_H3K27me3_3_CUTnRUN.filt.nodup.bam" #define path to my input bam file
new_input_file=${input_file%.bam}
sorted_bam=${input_file%.bam}.sorted.bam

#sort BAM file by name (sort -n), necessary for downstream analysis
#this step is necessary because we start from nodup.filter.BAM file took from Chip-seq piepline and it is sorted by position
#we need to sort it again by name
sambamba sort -n -t 4 -o $sorted_bam $input_file

#modification BAM to BED file (bedpe = paired end)
bedtools bamtobed -bedpe -i $sorted_bam > ${new_input_file}.bed

#filtering reads: keep only those mapping in same chr ($1==$4) and fragment distance <1000 bases ($6-$2 < 1000)
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${new_input_file}.bed > ${new_input_file}.clean.bed

#take only these fields: chr, start from one read, end from the paired read
#sort it by chr, by start position and by end position (in that order)
cut -f 1,2,6 ${new_input_file}.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${new_input_file}.fragments.bed

#transform to bedgraph (-bg) using a reference genome (-g)
bedtools genomecov -bg -i ${new_input_file}.fragments.bed -g /mnt/beegfs/rdeharo/prueba/Homo_sapiens.GRCh38.104.dna.all_chr.fa.fai > ${new_input_file}.fragments.bedgraph
```
