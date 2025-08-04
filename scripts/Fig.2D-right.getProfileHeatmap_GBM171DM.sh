#!/bin/bash
export pref='B171.HD2.DiffBind'
export bwdir='/./data/misc/GBM171DM/chipseq/'
export bwsuf='all_treat_pileup.bw'

run1Mock=A0088
run1HD60=A0090
run1HD120=A0093

export mock1=${bwdir}${run1Mock}.${bwsuf}
export hd_2a=${bwdir}${run1HD60}.${bwsuf}
export hd_3a=${bwdir}${run1HD120}.${bwsuf}

#cluster job setting parameters:
NTHREAD=8 #cpu usage per job

##########################################################################################
maindir=$PWD
regionDir=/./data/misc/GBM171DM/chipseq/

export extbp=1200
export deeptools='singularity run ./singularity/deeptools_3.4.3.sif'
export calcMatrix="$deeptools computeMatrix reference-point --referencePoint center -p $NTHREAD --beforeRegionStartLength $extbp --afterRegionStartLength $extbp --missingDataAsZero --skipZeros -R $regionDir"

#------------------------------------

function odir {
outputDir=$1
if [ ! -d $outputDir ]; then
    echo "Create $outputDir"
    mkdir $outputDir
fi
}

function pbsjob {
    RUN=$1
    FDR=$2

jobname=${RUN}.${FDR}.getHeatmap
pbs="$jobname.sge"
outMatfile=$pref.${RUN}.${FDR}.signal.${extbp}bp.gz

cat << ENDHERE > $pbs

#!/bin/bash
#$ -q analysis.q
#$ -cwd
#$ -l mfree=10G
#$ -l d_rt=4:0:0
#$ -S /bin/bash
#$ -N $jobname
#$ -o $jobname.sge.o
#$ -e $jobname.sge.e

##################

${calcMatrix}/$pref.${RUN}.${FDR}.bed -S $mock1 $hd_2a $hd_3a -o $outMatfile
$deeptools plotHeatmap -m $outMatfile -o $pref.${RUN}.${FDR}.signal.${extbp}bp.heatmap.pdf --samplesLabel $run1Mock  $run1HD60 $run1HD120 --whatToShow 'heatmap and colorbar' --regionsLabel "" --colorMap RdYlBu_r --zMax 8

date
ENDHERE

qsub $pbs
cd $maindir 

}

#------- loop & submit
for reg in chr ; do
    for fdr in FDRth0.01 ;  do
       pbsjob  $reg $fdr ;
    done
    sleep 1
done
