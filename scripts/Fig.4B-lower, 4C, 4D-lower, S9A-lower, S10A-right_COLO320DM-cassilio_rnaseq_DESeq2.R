
#input files
libNames = read.table('/net/nwgc/vol1/nobackup/nocleanup/tungch/RNA_AT4-48/rnaseq.fastqlist.txt', stringsAsFactors=F)[,2]
libNames <- libNames[c(31:33, 1:30,34:36)]
ucond = c('EV', 'SE1', 'SE2','SE3', 'SE4','SE7','SE9', 'Combo_SE1to4','Combo_CTL_SE7to9','Combined_SE1to4','Combined_CTL_SE7to9','untransfected')
conds = rep(ucond, each=3)
repli = rep(1:3, length(ucond))
names(conds) = libNames

rnaseqDir = "/net/nwgc/vol1/nobackup/nocleanup/tungch/RNA_AT4-48/"
tssbed='/net/nwgc/vol1/sharing/Wei_Lab/proj-in-situ-chia-pet/annotation/hg38.Genes.itxAnnotation.bed'
geneSet = '.gencode.v36.chialinPick'

sampleFiles = paste0(libNames, "/countGene/", libNames, geneSet, '.htseqcount')
sampleTable <- data.frame(sampleName = libNames,
                          fileName = sampleFiles,
                          condition = factor(conds, level=ucond),
                          rep = repli
                          )

sampleTable <- sampleTable[-c(16:21,25:27,31:36),]

########################################
# Load Libraries
########################################
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(dplyr)
library(data.table)
library(tidyr)
library(pheatmap)

###########################################################################################

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

shortnames = paste0(conds, '_',sub('RNA_AT4-48', '',libNames))
names(shortnames) = libNames

normCounts = counts(dds, normalized=TRUE)

normCounts.1 <- normCounts
colnames(normCounts.1) =  shortnames[colnames(normCounts.1)]
cormat = cor(normCounts.1) 
ph = pheatmap(cormat, 
		cluster_rows = F,
		cluster_cols = F)

library(ggplot2)

ggsave(ph,filename = "cor_AT4-48_ev_woCtrl_ordered.pdf")

PCA_vsdata <-  plotPCA(vsdata)

ggsave(PCA_vsdata,filename = "PCA_AT4-48_ev_woCtrl.pdf")

colnames(normCounts) = paste0(colnames(normCounts), "_norm")
normCounts.df = data.frame(geneName = geneDict[rownames(normCounts)], normCounts, check.names=F)

# comparisons

res1 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('condition', 'SE1', 'EV'))
res1.LFC <- lfcShrink(dds, coef="condition_SE1_vs_EV",, type="apeglm", res=res1)
res2 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('condition', 'SE2', 'EV'))
res2.LFC <- lfcShrink(dds, coef="condition_SE2_vs_EV", type="apeglm", res=res2)
res3 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('condition', 'SE3', 'EV'))
res3.LFC <- lfcShrink(dds, coef="condition_SE3_vs_EV", type="apeglm", res=res3)
res4 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('condition', 'SE4', 'EV'))
res4.LFC <- lfcShrink(dds, coef="condition_SE4_vs_EV", type="apeglm", res=res4)

res_Combo_SE1to4 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('condition', 'Combo_SE1to4', 'EV'))
res_Combo_SE1to4.LFC <- lfcShrink(dds, coef="condition_Combo_SE1to4_vs_EV", type="apeglm", res=res_Combo_SE1to4)
res_Combined_SE1to4 <- results(dds, alpha=0.05, lfcThreshold=logFCth, contrast=c('condition', 'Combined_SE1to4', 'EV'))
res_Combined_SE1to4.LFC <- lfcShrink(dds, coef="condition_Combined_SE1to4_vs_EV", type="apeglm", res=res_Combined_SE1to4)

#rename cols then add ENS in each tables 
colnames(res1.LFC)[2:5] = paste0(colnames(res1.LFC)[2:5],'_SE1')
colnames(res2.LFC)[2:5] = paste0(colnames(res2.LFC)[2:5],'_SE2')
colnames(res3.LFC)[2:5] = paste0(colnames(res3.LFC)[2:5],'_SE3')
colnames(res4.LFC)[2:5] = paste0(colnames(res4.LFC)[2:5],'_SE4')

colnames(res_Combo_SE1to4.LFC)[2:5] = paste0(colnames(res_Combo_SE1to4.LFC)[2:5],'_Combo_SE1to4')
colnames(res_Combined_SE1to4.LFC)[2:5] = paste0(colnames(res_Combined_SE1to4.LFC)[2:5],'_Combined_SE1to4')

normCounts.df$ENS  = rownames(normCounts.df) #after checking they are unique
#Check: identical(rownames(resCombo.LFC), normCounts.df$ENS)

j = 2:5 #selected columns
MatDiff = data.frame(normCounts.df, res1.LFC[,j], res2.LFC[,j], res3.LFC[,j], res4.LFC[,j],res_Combo_SE1to4.LFC[,j],res_Combined_SE1to4.LFC[,j], check.names=F)
projname = 'RNA_AT4-48'
write.table(MatDiff, file=paste0(projname, '.rnaSeq.DE_emptyVec_woCtrlSE.MasterTable.tsv'), sep='\t', quote=F, col.names=T, row.names=F)
write.table(sampleTable, file=paste0(projname, '.sampleInfo.txt'),  sep="\t", quote=F, col.names=T, row.names=F)

#plot selected  genes
col1="#8a7066"
col2="#6298f0"
col3="#dd62f0"
col4="#eb1337"
groupColors = c(col1, col2, col3, col4)

#####################################################
###  Volcano Plots
#####################################################
#fdrGene = 1e-6
thpadj = 0.05; lfc = 1
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


compname='SE3_v_EV'
pdedf3 = subset(MatDiff, padj_SE3 < thpadj)
summary(abs(pdedf3$log2FoldChange_SE3) )

n.df = nrow(pdedf3)
upth = ((pdedf3$log2FoldChange_SE3 ) > lfc & pdedf3$padj_SE3  < thpadj)
downth = ((pdedf3$log2FoldChange_SE3 ) < lfc & pdedf3$padj_SE3  < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedf3$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

volc3 = ggplot(data=pdedf3, aes(log2FoldChange_SE3, -log10(padj_SE3))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))


compname='SE4_v_EV'
pdedf4 = subset(MatDiff, padj_SE4 < thpadj)
summary(abs(pdedf4$log2FoldChange_SE4) )

n.df = nrow(pdedf4)
upth = ((pdedf4$log2FoldChange_SE4 ) > lfc & pdedf4$padj_SE4  < thpadj)
downth = ((pdedf4$log2FoldChange_SE4 ) < lfc & pdedf4$padj_SE4  < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedf4$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

volc4 = ggplot(data=pdedf4, aes(log2FoldChange_SE4, -log10(padj_SE4))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))


compname='Combo_SE1to4_v_EV'
pdedf_Combo_SE1to4 = subset(MatDiff, padj_Combo_SE1to4 < thpadj)
summary(abs(pdedf_Combo_SE1to4$log2FoldChange_Combo_SE1to4) )

n.df = nrow(pdedf_Combo_SE1to4)
upth = ((pdedf_Combo_SE1to4$log2FoldChange_Combo_SE1to4 ) > lfc & pdedf_Combo_SE1to4$padj_Combo_SE1to4  < thpadj)
downth = ((pdedf_Combo_SE1to4$log2FoldChange_Combo_SE1to4 ) < lfc & pdedf_Combo_SE1to4$padj_Combo_SE1to4  < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedf_Combo_SE1to4$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

volc_Combo_SE1to4 = ggplot(data=pdedf_Combo_SE1to4, aes(log2FoldChange_Combo_SE1to4, -log10(padj_Combo_SE1to4))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))


compname='Combined_SE1to4_v_EV'
pdedf_Combined_SE1to4 = subset(MatDiff, padj_Combined_SE1to4 < thpadj)
summary(abs(pdedf_Combined_SE1to4$log2FoldChange_Combined_SE1to4) )

n.df = nrow(pdedf_Combined_SE1to4)
upth = ((pdedf_Combined_SE1to4$log2FoldChange_Combined_SE1to4 ) > lfc & pdedf_Combined_SE1to4$padj_Combined_SE1to4  < thpadj)
downth = ((pdedf_Combined_SE1to4$log2FoldChange_Combined_SE1to4 ) < lfc & pdedf_Combined_SE1to4$padj_Combined_SE1to4  < thpadj)
iup = which(upth)
idown = which(downth)
th = rep('The rest', n.df)
th[iup] = paste0('Up:', length(iup))
th[idown] = paste0('Down:', length(idown))
pdedf_Combined_SE1to4$threshold = factor(th, level= c(paste0('Up:', length(iup)), paste0('Down:', length(idown)), 'The rest'))

volc_Combined_SE1to4 = ggplot(data=pdedf_Combined_SE1to4, aes(log2FoldChange_Combined_SE1to4, -log10(padj_Combined_SE1to4))) +
  geom_point(aes(col=threshold, shape=threshold)) +
  xlab("log2Fold") + ylab("-log10(FDR)") +
  scale_color_manual(values=c("#e942f5cc","#5badd9cc","#c0c0c077")) +
  ggtitle(paste(compname)) + theme_bw(16) + guides(color=guide_legend(title=paste0(2^lfc,"fold", ";FDR<",thpadj)))  +
  scale_shape_manual(values = c(17,15,1)) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))


pdf(file=paste0(projname, '.wrtEV.volcanos.pdf'), width=8, height=6)
print(volc1 + geom_text_repel(data=subset(pdedf1, padj_SE1<ptext & abs(log2FoldChange_SE1)>fctext ), aes(label=geneName)))
print(volc2 + geom_text_repel(data=subset(pdedf2, padj_SE2<ptext & abs(log2FoldChange_SE2)>fctext ), aes(label=geneName)))
print(volc3 + geom_text_repel(data=subset(pdedf3, padj_SE3<ptext & abs(log2FoldChange_SE3)>fctext ), aes(label=geneName)))
print(volc4 + geom_text_repel(data=subset(pdedf4, padj_SE4<ptext & abs(log2FoldChange_SE4)>fctext ), aes(label=geneName)))
print(volc_Combo_SE1to4 + geom_text_repel(data=subset(pdedf_Combined_SE1to4, padj_Combined_SE1to4<ptext & abs(log2FoldChange_Combined_SE1to4)>fctext ), aes(label=geneName)))
print(volc_Combined_SE1to4 + geom_text_repel(data=subset(pdedf_Combined_SE1to4, padj_Combined_SE1to4<ptext & abs(log2FoldChange_Combined_SE1to4)>fctext ), aes(label=geneName)))

dev.off()

#ggsave(volc1,filename = "volcano_AT4-48_SE1_vs_ev.pdf")
#ggsave(volc2,filename = "volcano_AT4-48_SE2_vs_ev.pdf")
#ggsave(volc3,filename = "volcano_AT4-48_SE3_vs_ev.pdf")
#ggsave(volc4,filename = "volcano_AT4-48_SE4_vs_ev.pdf")
#ggsave(volc7,filename = "volcano_AT4-48_SE7_vs_ev.pdf")
#ggsave(volc9,filename = "volcano_AT4-48_SE9_vs_ev.pdf")

#ggsave(volc_Combo_SE1to4,filename = "volcano_AT4-48_Combo_SE1to4_vs_ev.pdf")
#ggsave(volc_Combo_CTL_SE7to9,filename = "volcano_AT4-48_Combo_CTL_SE7to9_vs_ev.pdf")
#ggsave(volc_Combined_SE1to4,filename = "volcano_AT4-48_Combined_SE1to4_vs_ev.pdf")
#ggsave(volc_Combined_CTL_SE7to9,filename = "volcano_AT4-48_Combined_CTL_SE7to9_vs_ev.pdf")

#############
#*annotation
#############

#load(MatDiff.RData)
library(dplyr)

geneCoord = read.delim('/net/nwgc/vol1/sharing/Wei_Lab/USERS/tungch/refs/hg38/gencode.v36.chialinPick.geneCoord.txt', stringsAsFactors=F)
colnames(geneCoord)[4:5] = paste0('gene',colnames(geneCoord)[4:5])

cellname1='COLO320DM'
fin1 = paste0('/net/nwgc/vol1/sharing/Wei_Lab/USERS/tungch/tjongh/ec/combine2Reps/',cellname1, '.chiadropGeneList.tsv')

gc1 = read.delim(fin1, stringsAsFactors=F, check.names=F)

#MatDiff = read.delim('AT4-48.rnaSeq.DE_emptyVec_woCtrlSE.MasterTable.tsv', stringsAsFactors=F)

m = merge(MatDiff, geneCoord, by=c('ENS', 'geneName'), all.x = T)

#rearrange columns
jcoord = (ncol(m) - 4):ncol(m)
m = m[,c(1:2,jcoord,3:(jcoord[1]-1))]

jcd = c(2,8:ncol(gc1)) #index chiadrop columns
m1 = merge(m, gc1[,jcd], by='ENS', all.x = T)
write.table(m1, file=paste0(cellname1, '.AT4-48.rnaSeq.DE_emptyVec_woCtrlSE.MasterTable.annotated.tsv'), sep='\t', quote=F, col.names=T, row.names=F)

###########################
#Violin plot
##########################

library(scales) # map numeric value to color
library(RColorBrewer) #提供常见的配色方案

# reshape data in R
library(reshape2)
library(plyr)
library(ggplot2)
library(ggpubr)

cellname1='COLO320DM'
fin1 = paste0(cellname1, '.AT4-48.rnaSeq.DE_emptyVec_woCtrlSE.MasterTable.annotated.tsv')
m1 = read.delim(fin1, stringsAsFactors=F, check.names=F)

maxfdr=0.05
epsilon=1

colnames(m1)

b=subset(m1, padj_Combined_SE1to4 <= maxfdr, select=c(log2FoldChange_Combined_SE1to4, ec, nGEM))

b$ec_tag <- ifelse(b$ec == FALSE|is.na(b$ec) == TRUE, FALSE, TRUE) 
b$group <- ifelse(b$ec_tag == TRUE, "Genes_on_ecDNA",
                  ifelse(is.na(b$nGEM) == TRUE, "non_ecDNA_targets", "ecDNA_targets"))

summary(b)

my_comparisons <- list( c("Genes_on_ecDNA", "non_ecDNA_targets"),c("Genes_on_ecDNA", "ecDNA_targets"),c("ecDNA_targets", "non_ecDNA_targets"))

my_comparisons <- list( c("ecDNA_targets", "non_ecDNA_targets"))

FC=paste0(colnames(b)[1],"_vs_EV")
group=paste0(colnames(b)[5])

stat_box_data <- function(y, upper_limit = max(iris$Sepal.Length) * 1.2) {
  return(
    data.frame(
      y = 0.98 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'median =', round(median(y), 2), '\n')
    )
  )
}

#violin plot

violin <- ggplot(b, aes(x=factor(b[,5],level=c("Genes_on_ecDNA","ecDNA_targets", "non_ecDNA_targets")), y=b[,1],fill=b[,5])) +
  coord_cartesian(ylim = c(-8, 18)) +
  scale_y_continuous(breaks=c(-5,0,5, 10, 15)) +
  scale_fill_brewer(palette="Blues",name = group) +
  geom_violin(bw =0.6) +
  geom_boxplot(width=0.2, fill="white", size = 0.3,colour = "black") +
               #,outlier.shape = NA) +
  theme(legend.position="right") +
  labs(title=paste0(FC,"_FDR",maxfdr),x="ecDNA targets defined by ChIA-Drop", y = "Log2FC") +
  #scale_fill_discrete(name = "New Legend Title") +
  theme_bw() +
  theme(legend.position="right",
        panel.border=element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
        #legend.text= element_text(face="bold", color="black",family = "Arial", size=24),
        plot.title = element_text(hjust = 0.5,size=18,face="bold"),
        axis.text.x = element_text(face="bold", color="black", size=12),
        axis.text.y = element_text(face="bold",  color="black", size=10),
        axis.title.x = element_text(face="bold", color="black", size=18),
        axis.title.y = element_text(face="bold",color="black", size=18)) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(max(iris$Sepal.Length) * 1.6,max(iris$Sepal.Length) * 1.8,max(iris$Sepal.Length) * 1.95),
                     method = "wilcox.test",
                     #method = "t.test",
                     method.args = list(alternative = "less"),
                     #label = "p.signif")
                     #label = paste0("p = ", after_stat("p.format")))
                     label = aes(label = "p.format")) +
  stat_summary(fun.data = stat_box_data,
               geom = "text",
               hjust = 0.5,
               vjust = -1)

violin

ggsave(violin,filename = paste0("violin_",FC,"_FDR",maxfdr,".pdf"))














