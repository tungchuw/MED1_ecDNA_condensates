#!/bin/bash
#$ -q analysis.q
#$ -cwd
#$ -l mfree=20G
#$ -l d_rt=6:0:0
#$ -S /bin/bash
#$ -N COLO_2
#$ -o COLO_2.sge.o
#$ -e COLO_2.sge.e

##################
singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/singularity/deeptools_3.3.2.sif computeMatrix \
    reference-point --referencePoint center -R ACD0004_mergeFrags.bed -S ACD0004.R1.q30.lenBE50.MAJOR_CHROM.bw ACD0005.R1.q30.lenBE50.MAJOR_CHROM.bw A0026.all_treat_pileup.bw A0010.all_treat_pileup.bw\
    -p 16 --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --missingDataAsZero --skipZeros --blackListFileName hg38-blacklist.v2.bed \
    -o COLO.chiadrop_med1_H3K27ac_5000.gz

singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/singularity/deeptools_3.3.2.sif plotHeatmap -m COLO.chiadrop_med1_H3K27ac_5000.gz \
    -out COLO320DM.chiadrop_med1_H3K27ac_5000.heatmap.pdf --colorMap RdYlBu_r RdYlBu_r RdYlBu_r RdYlBu_r --samplesLabel "ACD0004" "ACD0005" "MED1" "H3K27Ac"\
    --whatToShow 'heatmap and colorbar' --regionsLabel "" \
    --zMax 22 18 0.5 0.6 --sortRegions descend --sortUsingSamples 1

date
