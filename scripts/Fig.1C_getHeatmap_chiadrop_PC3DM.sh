#!/bin/bash
#$ -q analysis.q
#$ -cwd
#$ -l mfree=20G
#$ -l d_rt=6:0:0
#$ -S /bin/bash
#$ -N PC3
#$ -o PC3.sge.o
#$ -e PC3.sge.e

##################
singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/singularity/deeptools_3.3.2.sif computeMatrix \
    reference-point --referencePoint center -R ACD0013_mergeFrags.bed -S ACD0013.R1.q30.lenBE50.MAJOR_CHROM.bw ACD0014.R1.q30.lenBE50.MAJOR_CHROM.bw A0038.all_treat_pileup.bw A0015.all_treat_pileup.bw \
    -p 16 --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --missingDataAsZero --skipZeros --blackListFileName hg38-blacklist.v2.bed \
    -o PC3.chiadrop_med1_H3K27ac_5000.gz

singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/singularity/deeptools_3.3.2.sif plotHeatmap -m PC3.chiadrop_med1_H3K27ac_5000.gz \
    -out PC3.chiadrop_med1_H3K27ac_5000.heatmap.pdf --colorMap RdYlBu_r RdYlBu_r RdYlBu_r RdYlBu_r --samplesLabel "ACD0013" "ACD0014" "MED1" "H3K27Ac" \
    --whatToShow 'heatmap and colorbar' --regionsLabel "" \
    --zMax 22 21 0.8 0.8 --sortRegions descend --sortUsingSamples 1

date
