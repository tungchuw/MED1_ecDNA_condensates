#Differential analysis for cut & tag
#use seurat4 environment for apeglm
# or without activating conda : Rseurat
library(dplyr)
library(GenomicRanges)
library(DESeq2)
library(ggplot2)
library(viridis)
library(stringr)

cellname1='PC3DM+'
expset = 'AT4-30'
#expset = 'AT4-31'
FDR= "0.05"
load(paste0(expset,'.fdr',FDR,'.countMatrix.rda'))

#commands:
perlannot='perl /net/nwgc/vol1/sharing/Wei_Lab/USERS/tungch/tjongh/ec/itxAnnotation.4peaks.pl'
bedtools =paste0('singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools')
itxfolder="/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/ChIA-PIPE/diffloop/"

#Files:
regiondir='/net/nwgc/vol1/nobackup/nocleanup/tungch/regions/'
ecDNAfile = paste0(regiondir, cellname1, '.ecDNA.bed')
SEfile  = paste0(regiondir,  cellname1, '.SErose.bed')
med1file = paste0(regiondir, cellname1, '.MED1overlap_peak.bed')
pol2file = paste0(regiondir, cellname1, '.Pol2.narrowPeak')
genebed='/net/nwgc/vol1/nobackup/nocleanup/tungch/proj-in-situ-chia-pet/annotation/hg38.Genes.itxAnnotation.bed'
anchorfile = paste0('/net/nwgc/vol1/sharing/Wei_Lab/USERS/tungch/tjongh/ec/combine2Reps/',cellname1, '.anchors_ecGEMs.annot.tsv')
cdrop = read.delim(anchorfile, stringsAsFactors=F, check.names=F) #will annotate final master table with nGEM

ncdrop = nrow(cdrop)
condition = factor(sampleInfo$V3, levels=unique(sampleInfo$V3))
selectR = which(rowSums(countMat) > nlib) ## remove low count
dataS = countMat[selectR,]
selectPeaks = masterPeak[selectR]

#DEseq
dds = DESeqDataSetFromMatrix(countData = dataS,
                             colData = DataFrame(condition),
                             design = ~ condition)
DDS = DESeq(dds)
#Check: resultsNames(DDS)

#saved image
#save.image(paste0(expset,'.fdr',FDR, '_withoutCtrl.DDS.RData'))
#load(paste0(expset,'.fdr',FDR, '.DDS.RData'))

#
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")


#correlation plot
cormat = cor(normDDS) 

library(pheatmap)
ph = pheatmap(cormat,
                cluster_rows = F,
                cluster_cols = F)
library(ggplot2)

ggsave(ph,filename = paste0("cor_",expset,'.fdr',FDR,"_ev_ordered.pdf"))



res1 <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'SE1', 'Empty'))
res1.LFC <- lfcShrink(DDS, coef="condition_SE1_vs_Empty", type="apeglm", res=res1)
res2 <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'SE2', 'Empty'))
res2.LFC <- lfcShrink(DDS, coef="condition_SE2_vs_Empty", type="apeglm", res=res2)
res6 <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'SE6', 'Empty'))
res6.LFC <- lfcShrink(DDS, coef="condition_SE6_vs_Empty", type="apeglm", res=res6)
res7 <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'SE7', 'Empty'))
res7.LFC <- lfcShrink(DDS, coef="condition_SE7_vs_Empty", type="apeglm", res=res7)
resCombo <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'Combo', 'Empty'))
resCombo.LFC <- lfcShrink(DDS, coef="condition_Combo_vs_Empty", type="apeglm", res=resCombo)
resUntr  <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'Untr', 'Empty'))
resUntr.LFC <- lfcShrink(DDS, coef="condition_Untr_vs_Empty", type="apeglm", res=resUntr )


#rename the columns for each comparison
colnames(res1.LFC) = paste0(colnames(res1.LFC), '.SE1_v_EV')
colnames(res2.LFC) = paste0(colnames(res2.LFC), '.SE2_v_EV')
colnames(res6.LFC) = paste0(colnames(res6.LFC), '.SE6_v_EV')
colnames(res7.LFC) = paste0(colnames(res7.LFC), '.SE7_v_EV')
colnames(resCombo.LFC) = paste0(colnames(resCombo.LFC), '.Combo_v_EV')
colnames(resUntr.LFC) = paste0(colnames(resUntr.LFC), '.Untr_v_EV')
j = 2:5 #selected columns
MatDiff = data.frame(CHR=as.vector(seqnames(selectPeaks)), Start=start(selectPeaks), End=end(selectPeaks),  normDDS, res1.LFC[,j], res2.LFC[,j], res6.LFC[,j], res7.LFC[,j],resCombo.LFC[,j], resUntr.LFC[,j], check.names=F)
nreg = nrow(MatDiff)

save.image(paste0(expset,'.fdr',FDR, '_withoutCtrl.MatDiff.RData'))
#load(paste0(expset,'.fdr',FDR, '_withoutCtrl.MatDiff.RData'))

#cdropLab = rep('.',nreg) #to label what chiadrop anchor on cut_tag fragment
#Annotation 
tmpnameA = paste0(cellname1, '.tmpA.bed')
tmpnameB = paste0(cellname1, '.tmpB.bed')
tmpout  = paste0(cellname1, '.tmp.out')
write.table(cbind(MatDiff[,1:3],1:nreg  ), sep="\t", quote=F, col.names=F, row.names=F, file=tmpnameA)
jdrop=c(2:5,1)
write.table(cdrop[,jdrop], sep="\t", quote=F, col.names=F, row.names=F, file=tmpnameB)
#check: colnames(cdrop)[c(2:8,10:11,1)] #
cmd = paste(bedtools, 'intersect -loj -a', tmpnameA, '-b', tmpnameB, '|', bedtools, 'merge -i stdin -c 4,8,9 -o first,sum,distinct >', tmpout)
system(cmd)
o = read.table(tmpout, stringsAsFactors=F);
colnames(o) = c(colnames(MatDiff)[1:3], 'regID', 'nGEM', 'chiadropAnchor')

outdf = bind_cols(MatDiff[as.numeric(o$regID),],o[,(ncol(o)-1):ncol(o)])
#outdf = bind_cols(MatDiff[as.numeric(o$regID),],o)
outdf$nGEM[outdf$nGEM=='.'] = 0

message("Region files being used: \n ", pol2file , "\n ", ecDNAfile, "\n ", med1file, "\n ", SEfile )

nc = ncol(outdf)

#Start getting new annotations
tmpname = paste0(cellname1, '.tmp.bed')

SElabel = rep(0, nreg)
write.table(cbind(outdf[,1:3],1:nreg), sep="\t", quote=F, col.names=F, row.names=F, file=tmpname)
cmd = paste(bedtools, 'intersect -wo -a', tmpname, '-b', "PC3DM+.SErose.bed", '|', bedtools, 'merge -i stdin -c 4,8 -o first,distinct_sort_num >', tmpout)

system(cmd)
o = read.table(tmpout, stringsAsFactors=F)
#o1 = o[,c(1:4,8,9)]
#o = o1
colnames(o) = c('CHR','Start','End', 'index','SE')
SElabel[o$index] = o$SE; #SE ranking at the ordered region

head(o)

table(SElabel)

#Add gene & MED1 annotations
write.table(cbind(outdf[,c(1:3,1:3)],0,1:nreg,0), sep="\t", quote=F, col.names=F, row.names=F, file=tmpname)
cmd = paste(perlannot, tmpname,  genebed, pol2file, ecDNAfile, SEfile, med1file,  cellname1 )
system(cmd)
oitx = read.table(paste0(cellname1, '.annotated_itx.txt'), header=T, sep="\t", comment.char='', check.names=FALSE, stringsAsFactors=FALSE);
tmpsufs = c('.bedpe.anchors.bed', '.bedpe.anchors.gene.bed', '.cis.sigf.anchors.peak1.intersect.bed', '.cis.sigf.anchors.peak2.intersect.bed','.cis.sigf.anchors.peak3.intersect.bed', '.cis.sigf.anchors.peak4.intersect.bed', '.annotated_itx.txt')
system(paste('rm',  tmpnameA, tmpnameB, tmpname, tmpout,  paste0(cellname1,tmpsufs, collapse=' ')))

finaldf  = data.frame(outdf, ec=oitx$ecDNA=='LR', pol2Peak=oitx$Peak=='LR', rankSE=SElabel, MED1=oitx$MED1=='LR', PGI=oitx$FinalCodeL, TSS=oitx$TSS_L, GeneBody=oitx$GeneL, check.names=F, stringsAsFactors=F )
write.table(finaldf, file=paste0(expset,'_',cellname1,'.fdr',FDR,'.cut_tag.annotated.deSeqMasterTable.tsv'), sep='\t', quote=F, col.names=T, row.names=F)






