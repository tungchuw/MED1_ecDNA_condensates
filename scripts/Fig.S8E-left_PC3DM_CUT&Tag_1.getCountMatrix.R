#Differential analysis for cut & tag
#Use Rbase

library(dplyr)
library(GenomicRanges)
library(chromVAR) #working with bam

expset = 'AT4-30'
maindir=paste0('/net/nwgc/vol1/nobackup/nocleanup/tungch/cut_tag/',expset,'/')
peakdir=paste0('/net/nwgc/vol1/nobackup/nocleanup/tungch/cut_tag/cleanPeaks/',expset,'/')

FDR = 0.001
#expset = 'AT4-31'
#peakfile_suf='_peaks.bed' #original gopeak
peakfile_suf=paste0('.FDR',FDR,'.noBL.narrowPeak')
bamfile_suf='.hg38.F2048.q30.bam'
sampleInfo = read.table(paste0('sampleInfo.', expset, '.txt'), stringsAsFactors=F)

#Create a master peak list
libnames = sampleInfo$V1
nlib = length(libnames)
message(nlib, ' samples will be in the matrix')

mPeak = GRanges()
for (h in libnames ){
    peakfile=paste0(peakdir, h, peakfile_suf)
    peakRes = read.table(peakfile, stringsAsFactors=F)
    message(peakfile, ' read.')
    mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
}
masterPeak = reduce(mPeak)

#get fragment counts for each library based on the masterPeak
imat = 1
countMat = matrix(NA, length(masterPeak), nlib)
for (h in libnames ){
    bamfile=paste0(maindir, h, '/align/',h, bamfile_suf)
    frag_counts = getCounts(bamfile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam") 
    countMat[,imat] = counts(frag_counts)[,1]
    imat = imat + 1
}
colnames(countMat) = libnames

save.image(paste0(expset,'.fdr',FDR,'.countMatrix.rda'))
