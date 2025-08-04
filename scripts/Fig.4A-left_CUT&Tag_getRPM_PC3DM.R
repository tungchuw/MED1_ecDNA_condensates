#Differential analysis for cut & tag
#Use Rbase

library(dplyr)
library(GenomicRanges)
library(chromVAR) #working with bam
library(Rsamtools)
library(ggplot2)
library(reshape2)
library(tidyverse)

cellname = "PC3"
maindir='/net/nwgc/vol1/nobackup/nocleanup/tungch/cut_tag/AT4-30/'
peakdir=paste0('/net/nwgc/vol1/nobackup/nocleanup/tungch/cut_tag/cleanPeaks/AT4-30/')
expset = 'AT4-30'

bamfile_suf='.hg38.F2048.q30.bam'
sampleInfo = read.table(paste0('sampleInfo.', expset, '.txt'), stringsAsFactors=F)
colnames(sampleInfo) = c("Library","Rep","condition")

#Create a master peak list
libnames = sampleInfo$V1
nlib = length(libnames)
message(nlib, ' samples will be in the matrix')
mPeak = GRanges()
regions = read.table("/net/nwgc/vol1/nobackup/nocleanup/tungch/test/cut_get_RPM_SEs/PC3/regions.bed", stringsAsFactors=F)
message(peakfile, ' read.')
mPeak = GRanges(seqnames = regions$V1, IRanges(start = regions$V2, end = regions$V3), strand = "*") %>% append(mPeak, .)
imat = 1
countMat = matrix(NA, length(mPeak), nlib)
totalReads <- numeric(nlib)  # 初始化存储总片段数的向量

for (h in libnames) {
    bamfile = paste0(maindir, h, '/align/', h, bamfile_suf)
    bam <- BamFile(bamfile)
    totalReads[imat] <- countBam(bam)$records
    frag_counts = getCounts(bamfile, mPeak, paired = TRUE, by_rg = FALSE, format = "bam") 
    countMat[, imat] = counts(frag_counts)[, 1]
    imat = imat + 1
}

# RPM matrix
rpmMat = sweep(countMat, 2, totalReads, "/") * 1e6

countMat = as.data.frame(countMat[,1:16])
rpmMat = as.data.frame(rpmMat[,1:16])

# Add rownames and colnames
rownames(countMat) <- regions$V4
colnames(countMat) <- libnames
rownames(rpmMat) <- regions$V4
colnames(rpmMat) <- libnames

# write table for rpm
rpmMat_t = as.data.frame(t(rpmMat))
rpmMat_t$Library = row.names(rpmMat_t)
rpmMat_t_master <- left_join(rpmMat_t, sampleInfo, by="Library")
colnames(rpmMat_t_master)
write.table(rpmMat_t_master[,c(5,7,6,1:4)]  , file=paste0(cellname,'.',expset, '.cut&tag.RPM.tsv'), sep='\t', quote=F, col.names=T, row.names=F)

# write table for rpm normalized to vector
new_rpmMat <- sweep(rpmMat, 1, rpmMat[, 1], "/")
new_rpmMat_t = as.data.frame(t(new_rpmMat))
new_rpmMat_t$Library = row.names(new_rpmMat_t)
new_rpmMat_t_master <- left_join(new_rpmMat_t, sampleInfo, by="Library")
colnames(new_rpmMat_t_master)
write.table(new_rpmMat_t_master[,c(5,7,6,1:4)]  , file=paste0(cellname,'.',expset, '.cut&tag.RPM.normalized.tsv'), sep='\t', quote=F, col.names=T, row.names=F)

