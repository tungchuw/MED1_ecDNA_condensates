#!/bin/bash
#$ -q analysis.q
#$ -cwd
#$ -l mfree=20G
#$ -l d_rt=6:0:0
#$ -S /bin/bash
#$ -N B171
#$ -o B171.sge.o
#$ -e B171.sge.e

##################
singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/singularity/deeptools_3.3.2.sif computeMatrix \
    reference-point --referencePoint center -R ACD0035_mergeFrags.bed -S ACD0034.R1.q30.lenBE50.MAJOR_CHROM.bw ACD0035.R1.q30.lenBE50.MAJOR_CHROM.bw A0041.all_treat_pileup.bw B171_H3K27ac_A0006.all_treat_pileup.bw\
    -p 16 --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --missingDataAsZero --skipZeros --blackListFileName hg38-blacklist.v2.bed \
    -o B171.chiadrop_med1_H3K27ac_5000.gz

singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/singularity/deeptools_3.3.2.sif plotHeatmap -m B171.chiadrop_med1_H3K27ac_5000.gz \
    -out B171.chiadrop_med1_H3K27ac_5000.heatmap.pdf --colorMap RdYlBu_r RdYlBu_r RdYlBu_r RdYlBu_r --samplesLabel "ACD0034" "ACD0035" "MED1" "H3K27Ac"\
    --whatToShow 'heatmap and colorbar' --regionsLabel "" \
    --zMax 9 10 0.9 0.7 --sortRegions descend --sortUsingSamples 1

date
