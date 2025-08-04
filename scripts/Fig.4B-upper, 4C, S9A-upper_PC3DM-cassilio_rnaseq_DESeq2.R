
#input files
libNames = read.table('../fastqR1.txt', stringsAsFactors=F)[,2]
libNames <- libNames[c(16:18, 1:15)]
ucond = c('EV', 'SE1', 'SE2','SE6','SE7', 'comboSE' )
conds = rep(ucond, each=3)
repli = rep(1:3, length(ucond))
names(conds) = libNames

rnaseqDir = "/projects/wei-lab/proj-SuperEnhancer/rnaseq/transfection/"
tssbed='/projects/wei-lab/proj-in-situ-chia-pet/annotation/hg38.Genes.itxAnnotation.bed'
geneSet = '.gencode.v36.chialinPick'

sampleFiles = paste0(libNames, "/countGene/", libNames, geneSet, '.htseqcount')
sampleTable <- data.frame(sampleName = libNames,
                          fileName = sampleFiles,
                          condition = factor(conds, level=ucond),
                          rep = repli
                          )

########################################
# Load Libraries
########################################
library(DESeq2)
library(ggplot2)
library(ggrepel)
library (reshape2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(dplyr)
library(data.table)
library(tidyr)
library(pheatmap)

###########################################################################################
# RNA-Seq 
tmp = read.delim(paste0(rnaseqDir, sampleTable$fileName[1]), header=F, stringsAsFactors=F)
iexcl = grep('^__', tmp$V1); tmp <- tmp[-iexcl,]
geneDict <- tmp$V2
names(geneDict) = tmp$V1

#######################################
# Make dds using DESeqDataSetFromHTSeqCount
#######################################
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
                sampleTable = sampleTable,
                directory = rnaseqDir,
                design= ~ condition)
ddsOrig <- DESeq(ddsHTSeq)
rawcount = counts(ddsOrig)
htseq =  as.data.frame(rawcount)

nlib = length(libNames)
thcount = nlib #Gene count total threshold across all samples
logFCth = 1 #log2FoldChange threshold
imin <- rowSums(counts(ddsOrig)) > thcount

dds <- ddsOrig[imin,] #We will use this one for downstream

#just checking the number of genes
nrow(ddsOrig)
nrow(dds)

vsdata = vst(dds, blind=FALSE)

shortnames = paste0(conds, '_',sub('AT4-31-', '',libNames))
names(shortnames) = libNames

normCounts = counts(dds, normalized=TRUE)

normCounts.1 <- normCounts
colnames(normCounts.1) =  shortnames[colnames(normCounts.1)]
cormat = cor(normCounts.1) 
ph = pheatmap(cormat, cluster_rows = FALSE, cluster_cols = FALSE)
colnames(normCounts) = paste0(colnames(normCounts), "_norm")
normCounts.df = data.frame(geneName = geneDict[rownames(normCounts)], normCounts, check.names=F)

# comparisons

res1 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('condition', 'SE1', 'EV'))
res1.LFC <- lfcShrink(dds, coef="condition_SE1_vs_EV", type="apeglm", res=res1)
res2 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('condition', 'SE2', 'EV'))
res2.LFC <- lfcShrink(dds, coef="condition_SE2_vs_EV", type="apeglm", res=res2)
res6 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('condition', 'SE6', 'EV'))
res6.LFC <- lfcShrink(dds, coef="condition_SE6_vs_EV", type="apeglm", res=res6)
res7 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('condition', 'SE7', 'EV'))
res7.LFC <- lfcShrink(dds, coef="condition_SE7_vs_EV", type="apeglm", res=res7)
resCombo <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('condition', 'comboSE', 'EV'))
resCombo.LFC <- lfcShrink(dds, coef="condition_comboSE_vs_EV", type="apeglm", res=resCombo)

#rename cols then add ENS in each tables 
colnames(res1.LFC)[2:5] = paste0(colnames(res1.LFC)[2:5],'_SE1')
colnames(res2.LFC)[2:5] = paste0(colnames(res2.LFC)[2:5],'_SE2')
colnames(res6.LFC)[2:5] = paste0(colnames(res6.LFC)[2:5],'_SE6')
colnames(res7.LFC)[2:5] = paste0(colnames(res7.LFC)[2:5],'_SE7')
colnames(resCombo.LFC)[2:5] = paste0(colnames(resCombo.LFC)[2:5],'_comboSE')

normCounts.df$ENS  = rownames(normCounts.df) #after checking they are unique
#Check: identical(rownames(resCombo.LFC), normCounts.df$ENS)

j = 2:5 #selected columns
MatDiff = data.frame(normCounts.df, res1.LFC[,j], res2.LFC[,j], res6.LFC[,j], res7.LFC[,j],resCombo.LFC[,j], check.names=F)
projname = 'AT4-31'
write.table(MatDiff, file=paste0(projname, '.rnaSeq.DE_emptyVec.MasterTable.tsv'), sep='\t', quote=F, col.names=T, row.names=F)
write.table(sampleTable, file=paste0(projname, '.sampleInfo.txt'),  sep="\t", quote=F, col.names=T, row.names=F)

#plot selected  genes
col1="#8a7066"
col2="#6298f0"
col3="#dd62f0"
col4="#eb1337"
groupColors = c(col1, col2, col3, col4)

pdf(file=paste0(projname, '.wrtEV.corHeatmap_ordered.pdf'), width=6, height=5.4)
print(ph)
print(plotPCA(vsdata))
dev.off()

#####################################################
###  Volcano Plots
#####################################################
#fdrGene = 1e-6
thpadj = 0.01; lfc = 1
ptext  = 1e-5; fctext = 4

compname='SE1_v_EV'
pdedf1 = subset(MatDiff, padj_SE1 < thpadj)
summary(abs(pdedf1$log2FoldChange_SE1) )

n.df = nrow(pdedf1)
upth = ((pdedf1$log2FoldChange_SE1 ) > lfc & pdedf1$padj_SE1  < thpadj)
downth = ((pdedf1$log2FoldChange_SE1 ) < lfc & pdedf1$padj_SE1  < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedf1$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

volc1 = ggplot(data=pdedf1, aes(log2FoldChange_SE1, -log10(padj_SE1))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))

compname='SE2_v_EV'
pdedf2 = subset(MatDiff, padj_SE2 < thpadj)
summary(abs(pdedf2$log2FoldChange_SE2) )

n.df = nrow(pdedf2)
upth = ((pdedf2$log2FoldChange_SE2 ) > lfc & pdedf2$padj_SE2  < thpadj)
downth = ((pdedf2$log2FoldChange_SE2 ) < lfc & pdedf2$padj_SE2  < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedf2$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

volc2 = ggplot(data=pdedf2, aes(log2FoldChange_SE2, -log10(padj_SE2))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))


compname='SE6_v_EV'
pdedf6 = subset(MatDiff, padj_SE6 < thpadj)
summary(abs(pdedf6$log2FoldChange_SE6) )

n.df = nrow(pdedf6)
upth = ((pdedf6$log2FoldChange_SE6 ) > lfc & pdedf6$padj_SE6  < thpadj)
downth = ((pdedf6$log2FoldChange_SE6 ) < lfc & pdedf6$padj_SE6  < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedf6$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

volc6 = ggplot(data=pdedf6, aes(log2FoldChange_SE6, -log10(padj_SE6))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))


compname='SE7_v_EV'
pdedf7 = subset(MatDiff, padj_SE7 < thpadj)
summary(abs(pdedf7$log2FoldChange_SE7) )

n.df = nrow(pdedf7)
upth = ((pdedf7$log2FoldChange_SE7 ) > lfc & pdedf7$padj_SE7  < thpadj)
downth = ((pdedf7$log2FoldChange_SE7 ) < lfc & pdedf7$padj_SE7  < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedf7$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

volc7 = ggplot(data=pdedf7, aes(log2FoldChange_SE7, -log10(padj_SE7))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))


compname='comboSE_v_EV'
pdedfCombo = subset(MatDiff, padj_comboSE < thpadj)
summary(abs(pdedfCombo$log2FoldChange_comboSE) )

n.df = nrow(pdedfCombo)
upth = ((pdedfCombo$log2FoldChange_comboSE ) > lfc & pdedfCombo$padj_comboSE  < thpadj)
downth = ((pdedfCombo$log2FoldChange_comboSE ) < lfc & pdedfCombo$padj_comboSE  < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedfCombo$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

volcCombo = ggplot(data=pdedfCombo, aes(log2FoldChange_comboSE, -log10(padj_comboSE))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))



pdf(file=paste0(projname, '.wrtEV.volcanos.pdf'), width=8, height=6)
print(volc1 + geom_text_repel(data=subset(pdedf1, padj_SE1<ptext & abs(log2FoldChange_SE1)>fctext ), aes(label=geneName)))
print(volc2 + geom_text_repel(data=subset(pdedf2, padj_SE2<ptext & abs(log2FoldChange_SE2)>fctext ), aes(label=geneName)))
print(volc6 + geom_text_repel(data=subset(pdedf6, padj_SE6<ptext & abs(log2FoldChange_SE6)>fctext ), aes(label=geneName)))
print(volc7 + geom_text_repel(data=subset(pdedf7, padj_SE7<ptext & abs(log2FoldChange_SE7)>fctext ), aes(label=geneName)))
print(volcCombo + geom_text_repel(data=subset(pdedfCombo, padj_comboSE<ptext & abs(log2FoldChange_comboSE)>fctext ), aes(label=geneName)))
dev.off()


