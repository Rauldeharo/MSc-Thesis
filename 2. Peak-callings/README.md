# Peak-callings used to obtain peak datasets

## MACS2
To obtain peaks using this peak-caller developed for ChIP-seq data, the most available version of MACS2 (MACS v2.2.8) was used.

All peak-datasets for each histone mark were obtained running MACS2 with an input control.

Peaks were called as following, using the Filtered BAM files sorted by position:

```
macs2 callpeak -t $BAM -c $BAM_INPUT -f $FORMAT -g $ORG -n ${PREFIX} --outdir $OUTDIR/peaks
```

## SEACR
Peak obtainment with this peak-caller was performed using the latest version (SEACR v1.3)

2 different modes were used (stringent and relaxed), which corresponds to a high or softer threshold in peak signal detection, respectively:

Using as input files: Experimental bedgraph ($exp_sample) and Control bedgraph ($IgG)
Obtaining as output $nB_WT_IgG_3_CUTnRUN

```
!/bin/bash

SEACR_1.3.sh $exp_sample $IgG norm stringent|relaxed nB_WT_IgG_3_CUTnRUN
```

