
#input files
runIDs = paste0('ARP00',44:51)
conds = rep(c('Mock', 'HD1h', 'HD2h','HD4h'), 2)
repli = rep(1:2, each=4)
tpoints = factor(rep(c(0,60,120,240), 2))
cellname='COLO320DM'
projname=paste0(cellname,'-HD2-RNAseq')
names(conds) = runIDs
cs = split(conds, conds)
cond.sorted = c(cs$Mock, cs$HD1h, cs$HD2h, cs$HD4h)

rnaseqDir = "/projects/wei-lab/proj-SuperEnhancer/rnaseq/HDtreatment/COLO320DM/"
tssbed='/projects/wei-lab/proj-in-situ-chia-pet/annotation/hg38.Genes.itxAnnotation.bed'
geneSet = '.gencode.v36.chialinPick'

sampleFiles = paste0(runIDs, "/countGene/", runIDs, geneSet, '.htseqcount')
sampleTable <- data.frame(sampleName = runIDs,
                          fileName = sampleFiles,
                          condition = factor(conds, level=conds[1:4]),
                          rep = repli,
                          timepoint = tpoints)

########################################
# Load Libraries
########################################
library(DESeq2)
library("ggplot2")
library("ggrepel")
library ("reshape2")
library(RColorBrewer)
library(gridExtra)
library(grid)
library(dplyr)
library(data.table)
library(tidyr)
library(pheatmap)

# ----------------- function -------------------

reorder_cormat <- function(cormat){
# Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd, method='ward.D')
    cormat <-cormat[hc$order, hc$order]
}

###########################################################################################

tmp = read.delim(paste0(rnaseqDir, sampleTable$fileName[1]), comment.char = "_", header=F, stringsAsFactors=F)
geneDict <- tmp$V2
names(geneDict) = tmp$V1

#######################################
# Make dds using DESeqDataSetFromHTSeqCount
#######################################
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
                sampleTable = sampleTable,
                directory = rnaseqDir,
                design= ~ timepoint)
ddsOrig <- DESeq(ddsHTSeq)
rawcount = counts(ddsOrig)
htseq =  as.data.frame(rawcount)
s = apply(htseq, 2, summary)

mypar=par
par(mfrow=c(4,2))
for (i in 1:length(cond.sorted)){
    l = names(cond.sorted)[i]
    hist(htseq[[l]], xlab='Raw count', breaks=400000, xlim=c(0,25), main=l)
}
par=mypar
sumcount = rowSums(htseq)
sumcount[sumcount > 1000000] = 1000000 #just for visualization at low range
sc1 = sumcount; sc1[sumcount < 1] = NA
sc2 = sumcount; sc2[sumcount < 2] = NA
sc3 = sumcount; sc3[sumcount < 3] = NA
sc4 = sumcount; sc4[sumcount < 4] = NA
sc5 = sumcount; sc5[sumcount < 5] = NA
sc6 = sumcount; sc6[sumcount < 6] = NA
sc10 = sumcount; sc10[sumcount < 10] = NA
par(mfrow=c(3,1))
hh = hist(sumcount, breaks=1000000, xlim=c(0,50), main='Histograms of sums over all 8 samples', xlab='Count per gene')
hist(sc1, breaks=1000000, xlim=c(0,50), main='Histograms of sums over all 8 samples', xlab='Count per gene')
hist(sc2, breaks=1000000, xlim=c(0,50), main='Histograms of sums over all 8 samples', xlab='Count per gene')
hist(sc3, breaks=1000000, xlim=c(0,50), main='Histograms of sums over all 8 samples', xlab='Count per gene')
hist(sc4, breaks=1000000, xlim=c(0,50), main='Histograms of sums over all 8 samples', xlab='Count per gene')
hist(sc5, breaks=1000000, xlim=c(0,50), main='Histograms of sums over all 8 samples', xlab='Count per gene')
hist(sc6, breaks=1000000, xlim=c(0,50), main='Histograms of sums over all 8 samples', xlab='Count per gene')
hist(sc10, breaks=1000000, xlim=c(0,50), main='Histograms of sums over all 8 samples', xlab='Count per gene')
hr = hh$counts[2:(ncount+1)]/hh$counts[1:ncount] #ratio of next count to the previous one

thcount = 1 #Gene count total threshold across all samples
logFCth = 1 #log2FoldChange threshold
imin <- rowSums(counts(ddsOrig)) > thcount
dds <- ddsOrig[imin,] #We will use this one for downstream

#just checking the number of genes
nrow(ddsOrig)
nrow(dds)

vsdata = vst(dds, blind=FALSE)
normCounts = counts(dds, normalized=TRUE)
cormat = cor(normCounts) 
cormat.ordered = cormat[c(rbind(1:4,1:4+4)),c(rbind(1:4,1:4+4))]
ph = pheatmap(cormat.ordered)
normCounts.df = data.frame(geneName = geneDict[rownames(normCounts)], normCounts)

# 1h vs 0 min comparison
res60v0 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('timepoint', '60', '0'))
res60v0LFC <- lfcShrink(dds, coef="timepoint_60_vs_0", type="apeglm", res=res60v0)

# 2 h vs 0 min comparison
res120v0<-results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('timepoint', '120', '0'))
res120v0LFC <- lfcShrink(dds, coef="timepoint_120_vs_0", type="apeglm", res=res120v0)
summary(res120v0LFC)

# 4 h vs 0 min comparison
res240v0<-results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('timepoint', '240', '0'))
res240v0LFC <- lfcShrink(dds, coef="timepoint_240_vs_0", type="apeglm", res=res240v0)
summary(res240v0LFC)

# check if gene response is influenced by timepoint
ddsLRT <- DESeq(dds, test="LRT", reduced=~1)
resLRT <- results(ddsLRT)

#rename cols then add ENS in each tables 
colnames(res60v0LFC)[2:5] = paste0(colnames(res60v0LFC)[2:5],'_60v0')
colnames(res120v0LFC)[2:5] = paste0(colnames(res120v0LFC)[2:5],'_120v0')
colnames(res240v0LFC)[2:5] = paste0(colnames(res240v0LFC)[2:5],'_240v0')
colnames(resLRT )[2:ncol(resLRT)] = paste0(colnames(resLRT)[2:ncol(resLRT)],'_LRT')
normCounts.df$Ens  = rownames(normCounts.df) #after checking they are unique
res60v0LFC$Ens  = rownames(res60v0LFC) 
res120v0LFC$Ens  = rownames(res120v0LFC) 
res240v0LFC$Ens  = rownames(res240v0LFC) 

#  write out the master table
res1 <- as.data.frame(res60v0LFC); res1$Ens <- rownames(res1);
res2 <- as.data.frame(res120v0LFC); res2$Ens <- rownames(res2);
res3 <- as.data.frame(res240v0LFC); res3$Ens <- rownames(res3);
res4 <- as.data.frame(resLRT); res3$Ens <- rownames(res3);
res.merged = merge(res1[,2:6], res2[,2:6], by='Ens')
res.merged = merge(res.merged, res3[,2:6], by='Ens')
outdf = merge(normCounts.df, res.merged, by='Ens')
outdf = outdf[,c(2:ncol(outdf),1)] 
write.table(sampleTable, file=paste0(projname, '.sampleInfo.txt'),  sep="\t", quote=F, col.names=T, row.names=F)
write.table(outdf, file=paste0(projname, '.rnaseqMasterTable.tsv'),  sep="\t", quote=F, col.names=T, row.names=F)

#plot selected  genes
col1="#8a7066"
col2="#6298f0"
col3="#dd62f0"
col4="#eb1337"
groupColors = c(col1, col2, col3, col4)

fdrGene = 1e-5
i = which(outdf$padj_60v0 < fdrGene);
sel60 = outdf[i[order(outdf$log2FoldChange_60v0[i])],]
sel60$geneName = factor(sel60$geneName, levels=sel60$geneName)
sel60.gathered = sel60[,1:(1+nrow(sampleTable))] %>% gather(colnames(sel60)[2:(1+nrow(sampleTable))], key='sample', value='NormalizedCounts')
sel60.gathered$Treatment=factor(conds[sel60.gathered$sample], level=conds[1:length(unique(conds))])

p60 = ggplot(sel60.gathered) + geom_point(aes(x=geneName, y=NormalizedCounts, color=Treatment)) + scale_y_log10() + xlab("Gene") + ylab("log10_NormCount") + ggtitle(paste0(nrow(sel60)," genes FDR_60v0<",fdrGene)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=groupColors)

i = which(outdf$padj_120v0 < fdrGene);
sel120 = outdf[i[order(outdf$log2FoldChange_120v0[i])],]
sel120$geneName = factor(sel120$geneName, levels=sel120$geneName)
sel120.gathered = sel120[,1:(1+nrow(sampleTable))] %>% gather(colnames(sel120)[2:(1+nrow(sampleTable))], key='sample', value='NormalizedCounts')
sel120.gathered$Treatment=factor(conds[sel120.gathered$sample], level=conds[1:length(unique(conds))])

p120 = ggplot(subset(sel120.gathered)) + geom_point(aes(x=geneName, y=NormalizedCounts, color=Treatment)) + scale_y_log10() + xlab("Gene") + ylab("log10_NormCount") + ggtitle(paste0(nrow(sel120)," genes FDR_120v0<",fdrGene)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=groupColors)

i = which(outdf$padj_240v0 < fdrGene);
sel240 = outdf[i[order(outdf$log2FoldChange_240v0[i])],]
sel240$geneName = factor(sel240$geneName, levels=sel240$geneName)
sel240.up = subset(sel240, log2FoldChange_240v0 > 0)
sel240.dn = subset(sel240, log2FoldChange_240v0 < 0)
sel240.gathered = sel240[,1:(1+nrow(sampleTable))] %>% gather(colnames(sel240)[2:(1+nrow(sampleTable))], key='sample', value='NormalizedCounts')
sel240.gathered$Treatment=factor(conds[sel240.gathered$sample], level=conds[1:length(unique(conds))])
sel240.gup = subset(sel240.gathered, geneName %in% sel240.up$geneName)
sel240.gdn = subset(sel240.gathered, geneName %in% sel240.dn$geneName)
p240 = ggplot(subset(sel240.gathered)) + geom_point(aes(x=geneName, y=NormalizedCounts, color=Treatment)) + scale_y_log10() + xlab("Gene") + ylab("log10_NormCount") + ggtitle(paste0(nrow(sel240)," genes FDR_240v0<",fdrGene)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=groupColors)
p240.up = ggplot(subset(sel240.gup)) + geom_point(aes(x=geneName, y=NormalizedCounts, color=Treatment)) + scale_y_log10() + xlab("Gene") + ylab("log10_NormCount") + ggtitle(paste0(nrow(sel240.up)," up genes FDR_240v0<",fdrGene)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=groupColors)
p240.dn = ggplot(subset(sel240.gdn)) + geom_point(aes(x=geneName, y=NormalizedCounts, color=Treatment)) + scale_y_log10() + xlab("Gene") + ylab("log10_NormCount") + ggtitle(paste0(nrow(sel240.dn)," down genes FDR_240v0<",fdrGene)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=groupColors)

pdf(file=paste0(projname, '.corHeatmap.pdf'), width=6, height=5.4)
print(ph)
print(plotPCA(vsdata, intgroup='timepoint'), ntop=1000)
dev.off()

pdf(file=paste0(projname, '.mostSignifGenes.pdf'), width=11, height=6)
print(p60)
print(p120)
print(p240)
print(p240.up)
print(p240.dn)
dev.off()

#####################################################
###  Volcano Plots
#####################################################
thpadj = 0.05; lfc = 1
ptext  = 1e-5; fctext = 2

pdedf = outdf
n.df = nrow(pdedf)
upth = ((pdedf$log2FoldChange_60v0) > lfc & pdedf$padj_60v0 < thpadj)
downth = ((pdedf$log2FoldChange_60v0) < lfc & pdedf$padj_60v0 < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedf$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

compname='60_vs_0'
volc = ggplot(data=pdedf, aes(log2FoldChange_60v0, -log10(padj_60v0))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))

#120_vs_0

pdedf120 = outdf
n.df = nrow(pdedf120)
upth = ((pdedf120$log2FoldChange_120v0) > lfc & pdedf120$padj_120v0 < thpadj)
downth = ((pdedf120$log2FoldChange_120v0) < lfc & pdedf120$padj_120v0 < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedf120$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

compname='120_vs_0'
volc120 = ggplot(data=pdedf120, aes(log2FoldChange_120v0, -log10(padj_120v0))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))


pdedf240 = outdf
n.df = nrow(pdedf240)
upth = ((pdedf240$log2FoldChange_240v0) > lfc & pdedf240$padj_240v0 < thpadj)
downth = ((pdedf240$log2FoldChange_240v0) < lfc & pdedf240$padj_240v0 < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedf240$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

compname='240_vs_0'
volc240 = ggplot(data=pdedf240, aes(log2FoldChange_240v0, -log10(padj_240v0))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))

pdf(file=paste0(projname, '.volcanos.pdf'), width=8, height=6)
print(volc + geom_text_repel(data=subset(pdedf, padj_60v0<ptext & abs(log2FoldChange_60v0)>fctext ), aes(label=geneName)))
print(volc120 + geom_text_repel(data=subset(pdedf120, padj_120v0<ptext & abs(log2FoldChange_120v0)>fctext ), aes(label=geneName)))
print(volc240 + geom_text_repel(data=subset(pdedf240, padj_240v0<ptext & abs(log2FoldChange_240v0)>fctext ), aes(label=geneName)))
dev.off()
