
#input files
runIDs = paste0('ARP00',52:57)
categories = c('Mock', 'HD2h', 'HD4h')
conds = rep(categories, each=2)
ncategory=length(categories)
repli = rep(1:2, ncategory)
tpoints = factor(rep(c(0,120,240), each=2))

cellname='B171'
projname=paste0(cellname,'-HD2-RNAseq')
names(conds) = runIDs
cs = split(conds, conds)
cond.sorted = c(cs$Mock, cs$HD2h, cs$HD4h)

rnaseqDir = paste0("/projects/wei-lab/proj-SuperEnhancer/rnaseq/HDtreatment/", cellname, '/')
tssbed='/projects/wei-lab/proj-in-situ-chia-pet/annotation/hg38.Genes.itxAnnotation.bed'
geneSet = '.gencode.v36.chialinPick'

sampleFiles = paste0(runIDs, "/countGene/", runIDs, geneSet, '.htseqcount')
sampleTable <- data.frame(sampleName = runIDs,
                          fileName = sampleFiles,
                          condition = factor(cond.sorted, levels=categories),
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
thcount = 1 #Gene count total threshold across all samples
logFCth = 1 #log2FoldChange threshold
imin <- rowSums(counts(ddsOrig)) > thcount
dds <- ddsOrig[imin,] #We will use this one for downstream
vsdata = vst(dds, blind=FALSE)
normCounts = counts(dds, normalized=TRUE)
cormat = cor(normCounts) 
ph = pheatmap(cormat, display_numbers = F, cluster_cols=F, cluster_rows=F)
normCounts.df = data.frame(geneName = geneDict[rownames(normCounts)], normCounts)

# 2 h vs 0 min comparison
res120v0<-results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('timepoint', '120', '0'))
res120v0LFC <- lfcShrink(dds, coef="timepoint_120_vs_0", type="apeglm", res=res120v0)

# 4 h vs 0 min comparison
res240v0<-results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('timepoint', '240', '0'))
res240v0LFC <- lfcShrink(dds, coef="timepoint_240_vs_0", type="apeglm", res=res240v0)

#rename cols then add ENS in each tables 
colnames(res120v0LFC)[2:5] = paste0(colnames(res120v0LFC)[2:5],'_120v0')
colnames(res240v0LFC)[2:5] = paste0(colnames(res240v0LFC)[2:5],'_240v0')
normCounts.df$Ens  = rownames(normCounts.df) #after checking they are unique
res120v0LFC$Ens  = rownames(res120v0LFC) 
res240v0LFC$Ens  = rownames(res240v0LFC) 

#  write out the master table
res2 <- as.data.frame(res120v0LFC); res2$Ens <- rownames(res2);
res3 <- as.data.frame(res240v0LFC); res3$Ens <- rownames(res3);
res.merged = merge(res2[,2:6], res3[,2:6], by='Ens')
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
fdrGene = 1e-2
i = which(outdf$padj_120v0 < fdrGene);
if(length(i) > 0){
sel120 = outdf[i[order(outdf$log2FoldChange_120v0[i])],]
sel120$geneName = factor(sel120$geneName, levels=sel120$geneName)
sel120.gathered = sel120[,1:(1+nrow(sampleTable))] %>% gather(colnames(sel120)[2:(1+nrow(sampleTable))], key='sample', value='NormalizedCounts')
sel120.gathered$Treatment=factor(conds[sel120.gathered$sample], level=categories)
p120 = ggplot(subset(sel120.gathered)) + geom_point(aes(x=geneName, y=NormalizedCounts, color=Treatment)) + scale_y_log10() + xlab("Gene") + ylab("log10_NormCount") + ggtitle(paste0(nrow(sel120)," genes FDR_120v0<",fdrGene)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=groupColors)

}

i = which(outdf$padj_240v0 < fdrGene);
if(length(i) > 0){
sel240 = outdf[i[order(outdf$log2FoldChange_240v0[i])],]
sel240$geneName = factor(sel240$geneName, levels=sel240$geneName)
sel240.up = subset(sel240, log2FoldChange_240v0 > 0)
sel240.dn = subset(sel240, log2FoldChange_240v0 < 0)
sel240.gathered = sel240[,1:(1+nrow(sampleTable))] %>% gather(colnames(sel240)[2:(1+nrow(sampleTable))], key='sample', value='NormalizedCounts')
sel240.gathered$Treatment=factor(conds[sel240.gathered$sample], level=categories)
sel240.gup = subset(sel240.gathered, geneName %in% sel240.up$geneName)
sel240.gdn = subset(sel240.gathered, geneName %in% sel240.dn$geneName)

p240 = ggplot(subset(sel240.gathered)) + geom_point(aes(x=geneName, y=NormalizedCounts, color=Treatment)) + scale_y_log10() + xlab("Gene") + ylab("log10_NormCount") + ggtitle(paste0(nrow(sel240)," genes FDR_240v0<",fdrGene)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=groupColors)
p240.up = ggplot(subset(sel240.gup)) + geom_point(aes(x=geneName, y=NormalizedCounts, color=Treatment)) + scale_y_log10() + xlab("Gene") + ylab("log10_NormCount") + ggtitle(paste0(nrow(sel240.up)," up genes FDR_240v0<",fdrGene)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=groupColors)
p240.dn = ggplot(subset(sel240.gdn)) + geom_point(aes(x=geneName, y=NormalizedCounts, color=Treatment)) + scale_y_log10() + xlab("Gene") + ylab("log10_NormCount") + ggtitle(paste0(nrow(sel240.dn)," down genes FDR_240v0<",fdrGene)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=groupColors)
}

pdf(file=paste0(projname, '.corHeatmap.pdf'), width=6, height=5.4)
print(ph)
print(plotPCA(vsdata, intgroup='timepoint'), ntop=1000)
dev.off()

pdf(file=paste0(projname, '.mostSignifGenes.pdf'), width=11, height=6)
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

##120_vs_0

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
#
compname='120_vs_0'
volc120 = ggplot(data=pdedf120, aes(log2FoldChange_120v0, -log10(padj_120v0))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))
#
#
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
#
compname='240_vs_0'
volc240 = ggplot(data=pdedf240, aes(log2FoldChange_240v0, -log10(padj_240v0))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))
#
#
pdf(file=paste0(projname, '.volcanos.pdf'), width=8, height=6)
print(volc120 + geom_text_repel(data=subset(pdedf120, padj_120v0<ptext & abs(log2FoldChange_120v0)>fctext ), aes(label=geneName)))
print(volc240 + geom_text_repel(data=subset(pdedf240, padj_240v0<ptext & abs(log2FoldChange_240v0)>fctext ), aes(label=geneName)))
dev.off()
