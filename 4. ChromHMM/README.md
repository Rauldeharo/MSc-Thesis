# ChromHMM
## Find out biological functionality on the peak-datasets
In order to find out biological functionality on these peak-datasets, a ChromHMM software (v1.24) was used. (http://compbio.mit.edu/ChromHMM/)

To run this software,  2 different commands were used for the following purposes:

### BinarizedBED | BinarizedBAM

This command was used to convert BED or BAM files into a binary files.

In this speicifc case, it is shown a ChromHMM for specific condition (ChromHMM directly from CUT&RUN BAM files)

```
#!/bin/bash

#we need to have all BAM or BED files in the same folder to indicate the directory

samples_directory="/mnt/beegfs/rdeharo/cutnrun_link_samples"
control_directory="/mnt/beegfs/rdeharo/cutnrun_link_controls"
chr_length_file="/mnt/beegfs/rdeharo/prueba/Homo_sapiens.GRCh38.104.dna.all_chr_2cols.fa.fai"
cell_mark_file="/mnt/beegfs/rdeharo/prueba/cell_mark_file_CnR.txt" #file indicating histone marks and cell type used
output_binary_directory="/mnt/beegfs/rdeharo/binarized_bam_CnR"

output_model_directory="/mnt/beegfs/rdeharo/ChromHMM_model_ChIP" #Directory to store ChromHMM Learned model

num_states=10 #Indicate number of states used for Learn model command

#Command to binarize BAM files
java -mx4000M -jar /opt/ohpc/pub/apps/ChromHMM/1.24/ChromHMM.jar BinarizeBam -paired -c $control_directory $chr_length_file $samples_directory $cell_mark_file $output_binary_directory

#File name modification for further analysis (binarized files needs to be named as cell_type_"chr"_number_"binary.txt"
cd $output_binary_directory
for file in $(ls -v $output_binary_directory);do
  number="${file#nB_WT_}"
  number="${number%_binary.txt}"
  sed -i "1s/.*/nB_WT\tchr$number/" $file
  mv $file "nB_WT_chr${number}_binary.txt"

done
```
### LearnModel

From previously binarized files, LearnModel command was used indicating the specific number of states to be learned in the model:

```
#Command to learn from binarized files and produce different chromatin states
java -mx20000M -jar /opt/ohpc/pub/apps/ChromHMM/1.24/ChromHMM.jar LearnModel -p 10 $output_binary_directory $output_model_directory $num_states hg38
```
