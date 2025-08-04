#Differential analysis for cut & tag
library(dplyr)
library(GenomicRanges)
library(DESeq2)
library(ggplot2)
library(viridis)
library(stringr)

cellname1='COLO320DM'
expset = 'C-n-T_AT4-48'
FDR= "0.01"
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

#saved image
save.image(paste0(expset,'.fdr',FDR, '_withoutCtrl.DDS.RData'))

#
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")

res1 <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'SE1', 'Empty'))
res1.LFC <- lfcShrink(DDS, coef="condition_SE1_vs_Empty", type="apeglm", res=res1)
res2 <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'SE2', 'Empty'))
res2.LFC <- lfcShrink(DDS, coef="condition_SE2_vs_Empty", type="apeglm", res=res2)
res3 <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'SE3', 'Empty'))
res3.LFC <- lfcShrink(DDS, coef="condition_SE3_vs_Empty", type="apeglm", res=res3)
res4 <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'SE4', 'Empty'))
res4.LFC <- lfcShrink(DDS, coef="condition_SE4_vs_Empty", type="apeglm", res=res4)

res_Combo_SE1.4 <- results(DDS, 
                        alpha=0.05, 
                        independentFiltering = FALSE, 
                        altHypothesis = "greaterAbs", 
                        contrast=c('condition', 'Combo_SE1.4', 'Empty'))
res_Combo_SE1.4.LFC <- lfcShrink(DDS, coef="condition_Combo_SE1.4_vs_Empty", type="apeglm", res=res_Combo_SE1.4)
res_Combined_SE1.4 <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'Combined_SE1.4', 'Empty'))
res_Combined_SE1.4.LFC <- lfcShrink(DDS, coef="condition_Combined_SE1.4_vs_Empty", type="apeglm", res=res_Combined_SE1.4)

resUntr  <- results(DDS, alpha=0.05, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c('condition', 'Untransfected', 'Empty'))
resUntr.LFC <- lfcShrink(DDS, coef="condition_Untransfected_vs_Empty", type="apeglm", res=resUntr )

#rename the columns for each comparison
colnames(res1.LFC) = paste0(colnames(res1.LFC), '.SE1_v_EV')
colnames(res2.LFC) = paste0(colnames(res2.LFC), '.SE2_v_EV')
colnames(res3.LFC) = paste0(colnames(res3.LFC), '.SE3_v_EV')
colnames(res4.LFC) = paste0(colnames(res4.LFC), '.SE4_v_EV')

colnames(res_Combo_SE1.4.LFC) = paste0(colnames(res_Combo_SE1.4.LFC), '.Combo_SE1.4_v_EV')
colnames(res_Combined_SE1.4.LFC) = paste0(colnames(res_Combined_SE1.4.LFC), '.Combined_SE1.4_v_EV')

colnames(resUntr.LFC) = paste0(colnames(resUntr.LFC), '.Untransfected_vs_Empty')

j = 2:5 #selected columns
MatDiff = data.frame(CHR=as.vector(seqnames(selectPeaks)), Start=start(selectPeaks), End=end(selectPeaks),  
                    normDDS, res1.LFC[,j], res2.LFC[,j], res3.LFC[,j], res4.LFC[,j], 
                    res_Combo_SE1.4.LFC[,j], res_Combined_SE1.4.LFC[,j], resUntr.LFC[,j], check.names=F)

colnames(MatDiff)

nreg = nrow(MatDiff)

save.image(paste0(expset,'.fdr',FDR, '_withoutCtrl.MatDiff.RData'))

#Annotation 
tmpnameA = paste0(cellname1, '.tmpA.bed')
tmpnameB = paste0(cellname1, '.tmpB.bed')
tmpout  = paste0(cellname1, '.tmp.out')
write.table(cbind(MatDiff[,1:3],1:nreg  ), sep="\t", quote=F, col.names=F, row.names=F, file=tmpnameA)
jdrop=c(2:5,1)
write.table(cdrop[,jdrop], sep="\t", quote=F, col.names=F, row.names=F, file=tmpnameB)
cmd = paste(bedtools, 'intersect -loj -a', tmpnameA, '-b', tmpnameB, '|', bedtools, 'merge -i stdin -c 4,8,9 -o first,sum,distinct >', tmpout)
system(cmd)
o = read.table(tmpout, stringsAsFactors=F);
colnames(o) = c(colnames(MatDiff)[1:3], 'regID', 'nGEM', 'chiadropAnchor')

outdf = bind_cols(MatDiff[as.numeric(o$regID),],o[,(ncol(o)-1):ncol(o)])
outdf$nGEM[outdf$nGEM=='.'] = 0
message("Region files being used: \n ", pol2file , "\n ", ecDNAfile, "\n ", med1file, "\n ", SEfile )
nc = ncol(outdf)

#Start getting new annotations
tmpname = paste0(cellname1, '.tmp.bed')
SElabel = rep(0, nreg)
write.table(cbind(outdf[,1:3],1:nreg), sep="\t", quote=F, col.names=F, row.names=F, file=tmpname)
cmd = paste(bedtools, 'intersect -wo -a', tmpname, '-b', "COLO320DM.SErose.bed", '|', bedtools, 'merge -i stdin -c 4,8 -o first,distinct_sort_num >', tmpout)
system(cmd)
o = read.table(tmpout, stringsAsFactors=F)
colnames(o) = c('CHR','Start','End', 'index','SE')
SElabel[o$index] = o$SE; #SE ranking at the ordered region

#Add gene & MED1 annotations
write.table(cbind(outdf[,c(1:3,1:3)],0,1:nreg,0), sep="\t", quote=F, col.names=F, row.names=F, file=tmpname)
cmd = paste(perlannot, tmpname,  genebed, pol2file, ecDNAfile, SEfile, med1file,  cellname1 )
system(cmd)
oitx = read.table(paste0(cellname1, '.annotated_itx.txt'), header=T, sep="\t", comment.char='', check.names=FALSE, stringsAsFactors=FALSE);
tmpsufs = c('.bedpe.anchors.bed', '.bedpe.anchors.gene.bed', '.cis.sigf.anchors.peak1.intersect.bed', '.cis.sigf.anchors.peak2.intersect.bed','.cis.sigf.anchors.peak3.intersect.bed', '.cis.sigf.anchors.peak4.intersect.bed', '.annotated_itx.txt')
system(paste('rm',  tmpnameA, tmpnameB, tmpname, tmpout,  paste0(cellname1,tmpsufs, collapse=' ')))

finaldf  = data.frame(outdf, ec=oitx$ecDNA=='LR', pol2Peak=oitx$Peak=='LR', rankSE=SElabel, MED1=oitx$MED1=='LR', PGI=oitx$FinalCodeL, TSS=oitx$TSS_L, GeneBody=oitx$GeneL, check.names=F, stringsAsFactors=F )
write.table(finaldf, file=paste0(expset,'_',cellname1,'.fdr',FDR,'.cut_tag.annotated.deSeqMasterTable_without_ctrlSE.tsv'), sep='\t', quote=F, col.names=T, row.names=F)
