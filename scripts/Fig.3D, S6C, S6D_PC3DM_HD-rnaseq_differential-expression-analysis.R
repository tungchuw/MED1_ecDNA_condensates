
#input files
runIDs = paste0('ARP00',38:43)
conds = rep(c('Mock', 'HD1h', 'HD2h'), 2)
reps = rep(1:2, each=3)
tpoints = factor(rep(c(0,60,120), 2))

projname='PC3DM+_HD2_RNAseq'
names(conds) = runIDs

rnaseqDir = "/projects/wei-lab/proj-SuperEnhancer/rnaseq/HDtreatment/PC3DM+/"
tssbed='./hg38.Genes.itxAnnotation.bed'
geneSet = '.gencode.v36.chialinPick'

sampleFiles = paste0(runIDs, "/countGene/", runIDs, geneSet, '.htseqcount')
sampleTable <- data.frame(sampleName = runIDs,
                          fileName = sampleFiles,
                          condition = factor(conds, level=conds[1:3]),
                          rep = reps,
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
# RNA-Seq 

tmp = read.delim(paste0(rnaseqDir, sampleFiles[1]), header=F, stringsAsFactors=F)
tmp = tmp[-grep('^__', tmp$V1),]
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
ph = pheatmap(cormat)
normCounts.df = data.frame(geneName = geneDict[rownames(normCounts)], normCounts)

# 1h vs 0 min comparison
res60v0 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('timepoint', '60', '0'))
res60v0LFC <- lfcShrink(dds, coef="timepoint_60_vs_0", type="apeglm", res=res60v0)
summary(res60v0LFC)

# 2 h vs 0 min comparison
res120v0<-results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('timepoint', '120', '0'))
res120v0LFC <- lfcShrink(dds, coef="timepoint_120_vs_0", type="apeglm", res=res120v0)

# check if gene response is influenced by timepoint
ddsLRT <- DESeq(dds, test="LRT", reduced=~1)
resLRT <- results(ddsLRT)

#rename cols then add ENS in each tables 
colnames(res60v0LFC)[2:5] = paste0(colnames(res60v0LFC)[2:5],'_60v0')
colnames(res120v0LFC)[2:5] = paste0(colnames(res120v0LFC)[2:5],'_120v0')
colnames(resLRT )[2:ncol(resLRT)] = paste0(colnames(resLRT)[2:ncol(resLRT)],'_LRT')
normCounts.df$Ens  = rownames(normCounts.df) #after checking they are unique
res60v0LFC$Ens  = rownames(res60v0LFC) 
res120v0LFC$Ens  = rownames(res120v0LFC) 

#  write out the master table
res1 <- as.data.frame(res60v0LFC); res1$Ens <- rownames(res1);
res2 <- as.data.frame(res120v0LFC); res2$Ens <- rownames(res2);
res3 <- as.data.frame(resLRT); res3$Ens <- rownames(res3);
res.merged = merge(res1[,2:6], res2[,2:6], by='Ens')

tmp = merge(normCounts.df, res.merged, by='Ens')
normCountsMetaTestTable = merge(tmp, res3, by='Ens')
outdf = normCountsMetaTestTable[,-c(18:22)]
write.table(sampleTable, file='trio.HDtreatment.txt',  sep="\t", quote=F, col.names=T, row.names=F)
outdf = outdf[,c(2:ncol(outdf), 1)]
write.table(outdf, file='PC3DM+_HD2.rnaseqMasterTable.tsv',  sep="\t", quote=F, col.names=T, row.names=F)

#plot selected  genes
col1="#8a7066"
col2="#6298f0"
col3="#dd62f0"

i = which(outdf$padj_60v0 < 1e-5);
sel60 = outdf[i[order(outdf$log2FoldChange_60v0[i])],]
sel60$geneName = factor(sel60$geneName, levels=sel60$geneName)
sel60.gathered = sel60[,1:7] %>% gather(colnames(sel60)[2:7], key='sample', value='NormalizedCounts')
sel60.gathered$Treatment=factor(conds[sel60.gathered$sample], level=conds[1:3])

p60 = ggplot(subset(sel60.gathered)) + geom_point(aes(x=geneName, y=NormalizedCounts, color=Treatment)) + scale_y_log10() + xlab("Gene") + ylab("log10_NormCount") + ggtitle(paste0(nrow(sel60)," genes FDR_60v0<1e-5")) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=c(col1, col2, col3))

i = which(outdf$padj_120v0 < 1e-5);
sel120 = outdf[i[order(outdf$log2FoldChange_120v0[i])],]
sel120$geneName = factor(sel120$geneName, levels=sel120$geneName)
sel120.gathered = sel120[,1:7] %>% gather(colnames(sel120)[2:7], key='sample', value='NormalizedCounts')
sel120.gathered$Treatment=factor(conds[sel120.gathered$sample], level=conds[1:3])

p120 = ggplot(subset(sel120.gathered)) + geom_point(aes(x=geneName, y=NormalizedCounts, color=Treatment)) + scale_y_log10() + xlab("Gene") + ylab("log10_NormCount") + ggtitle(paste0(nrow(sel120)," genes FDR_120v0<1e-5")) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=c(col1, col2, col3))

pdf(file=paste0(projname, '.corHeatmap.pdf'), width=6, height=5.5)
print(ph)
print(plotPCA(vsdata, intgroup='timepoint'), ntop=1000)
dev.off()

pdf(file=paste0(projname, '.mostSignifGenes.pdf'), width=11, height=6)
print(p60)
print(p120)
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


pdf(file=paste0(projname, '.volcanos.pdf'), width=8, height=6)
print(volc + geom_text_repel(data=subset(pdedf, padj_60v0<ptext & abs(log2FoldChange_60v0)>fctext ), aes(label=geneName)))
print(volc120 + geom_text_repel(data=subset(pdedf120, padj_120v0<ptext & abs(log2FoldChange_120v0)>fctext ), aes(label=geneName)))
dev.off()
