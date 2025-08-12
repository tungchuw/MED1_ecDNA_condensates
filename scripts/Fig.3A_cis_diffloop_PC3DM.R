
require(data.table)
require(ggplot2)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggrepel)

cellname='PC3DM+'
itxDir = "./ChIA-PIPE/diffloop/" #modify accordingly
samplePhenotype = read.table(paste0(itxDir,cellname, "_itxlist.txt")); 
colnames(samplePhenotype) = c("name","condition", "rep") #modify
celltreat <- sub('mock', '',paste(unique(samplePhenotype$condition), collapse='')) #mock must be in samplePhenotype
cellcntrl <- 'mock'
compname <- paste0(celltreat, '_vs_', cellcntrl)
nsample = nrow(samplePhenotype)

bedtools='singularity run ./singularity/chipseqtools.sif bedtools intersect -wo -a '
bedmerge='singularity run ./singularity/chipseqtools.sif bedtools merge -i '
bedpetools='singularity run ./singularity/chipseqtools.sif bedtools pairtopair -is -type both -a '

itxlibs = samplePhenotype$name;
itxsuf='.DA.itx'

#For annotation

perlannot='/usr/bin/perl ./itxAnnotation.4peaks.pl'
genebed='./hg38.Genes.itxAnnotation.bed'
regiondir='./data/misc/PC3DM/'

ecDNAfile=paste0(regiondir, cellname, '.ecDNA.bed')
peakfile=paste0(regiondir, cellname, ".Pol2.narrowPeak")
roseBed = paste0(regiondir, cellname , '.SErose.bed')
med1peak = paste0(regiondir, cellname, '.MED1overlap_peak.bed')

#-----------------------------------------------------------------------------------------------------------
#function
getAnchorIndex <- function(itx1, libID){
    n = nrow(itx1)
    anchorL = numeric(n);
    anchorR = numeric(n);
    names(anchorL) = names(anchorR) = 1:n
    tmpbed = 'tmp.bed'; tmpout = 'tmp.txt'

    m1 = itx1[,1:3]; m1$idL = 1:n
    write.table(m1, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
    com = paste0(bedtools, tmpbed, ' -b ', anchorfile, ' > ', tmpout)
    system(com)
    out = fread(tmpout, stringsAsFactors=F)
    #out = read.table(tmpout, stringsAsFactors=F)
    anchorL[out$V4] = out$V8;

    m1 = itx1[,4:6]; m1$idR = 1:n
    write.table(m1, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
    com = paste0(bedtools, tmpbed, ' -b ', anchorfile, ' > ', tmpout)
    system(com)
    out = fread(tmpout, stringsAsFactors=F)
    #out = read.table(tmpout, stringsAsFactors=F)
    anchorR[out$V4] = out$V8;
    rl = cbind(anchorL, anchorR)
    if(length(anchorL) != n){message('anchorL length not the same as: ',n)}
    if(length(anchorR) != n){message('anchorR length not the same as: ',n)}
    anchorL=apply(rl, 1, min)
    anchorR=apply(rl, 1, max)
    anchorLab = paste0( anchorL, '_', anchorR);
    idPETs = data.frame(pairID=anchorLab, iPET=itx1[,7], stringsAsFactors=FALSE)
    idPETs = aggregate(iPET~pairID, idPETs, sum)
    colnames(idPETs)[2] = paste0(libID, '_', colnames(idPETs)[2])
    return(idPETs)
}

#==============================
#Reading itx.DA
itxDir = "./ChIA-PIPE/HDtreat/" #modify accordingly
itxs = list()
for ( i in itxlibs ){
    f1 =  paste0(itxDir, i, "/", i, itxsuf)
    headerset= c('ChrL', 'StartL', 'EndL', 'ChrR', 'StartR', 'EndR', 'iPET')
    message('Read  ', f1)
    #m1 = read.table(f1, header=FALSE, stringsAsFactors=F)
    m1 = fread(f1, header=FALSE, stringsAsFactors=F)[,1:7]
    colnames(m1) = headerset
    itxs[[i]] = m1
}

allitx = bind_rows(itxs) #combine all itx
anchors1 = allitx[,1:3]
anchors2 = allitx[,4:6]
colnames(anchors1) = colnames(anchors2) = c('CHR','Start','End')
anchors = rbind(anchors1, anchors2)
rm(anchors1,anchors2)

#indexing the combined anchors
tmpbed = paste0(cellname, '.tmp.bed'); 
tmpout = paste0(cellname, '.tmp.out');
write.table(anchors, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0('sort -k1,1 -k2,2n ', tmpbed, ' | ', bedmerge, ' stdin > ', tmpout) 
system(com)
anchors = fread(tmpout, stringsAsFactors=F)
#colnames(anchors) =  c('CHR','Start','End'); anchors$ID = 1:nrow(anchors)
colnames(anchors) =  c('CHR','Start','End'); anchors$ID = paste0(anchors$CHR, ':', anchors$Start,'-',anchors$End)
anchorfile = paste0(cellname, '.anchors.bed')
write.table(anchors, file=anchorfile, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0(bedtools, anchorfile, ' -b ', ecDNAfile, ' > ', tmpout)
system(com)
o = read.table(tmpout, stringsAsFactors=F)

#annotate anchors for EC or not
atype = rep('ch', nrow(anchors))
atype[which(anchors$ID %in% o$V4)] = 'ec'
names(atype) = anchors$ID

anchors$locus <- atype;

#map the itx anchors with ID
anchor.list = list()
for ( i in itxlibs ){
    message('Getting pairID for ', i)
    anchor.list[[i]] = getAnchorIndex(itxs[[i]], libID=i)
}

#alltab = anchor.list %>% reduce(full_join, by="pairID") #merging the tables
a = merge(anchor.list[[1]], anchor.list[[2]], all=T)
b = merge(anchor.list[[3]], anchor.list[[4]], all=T)

alltab = merge(a,b, all=T)
alltab[is.na(alltab)] = 0;

#Writing input for DEseq2
for ( i in 1:4 ){
    lib1 = itxlibs[i]
    message('Writing labeled itx for ', lib1)
    write.table(alltab[,c(1,i+1)], file=paste0(lib1, '.petID.txt'), sep='\t', quote=F, col.names=F, row.names=F)
}



#DEseq approach
sampleFiles = paste0(samplePhenotype$name, '.petID.txt')
sampleTable <- data.frame(sampleName = samplePhenotype$name,
                          fileName = sampleFiles,
                          condition = factor(samplePhenotype$condition))

#######################################
# Make dds using DESeqDataSetFromHTSeqCount
#######################################
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       design= ~ condition)
# Specify a reference level
ddsHTSeq$condition = relevel(ddsHTSeq$condition, ref = 'mock')
ddsOrig = DESeq(ddsHTSeq)
#dds.de = DESeq(ddsHTSeq)
thcount <- 6
imin <- rowSums(counts(ddsOrig)) > thcount
dds.de <- ddsOrig[imin,] #We will use this one for downstream

summary(rowSums(counts(ddsOrig)))

#just checking the number of genes
nrow(ddsOrig)
nrow(dds.de)


#######################################
# Get normalized counts
#######################################
rawcounts =  counts(dds.de, normalized=F)
colnames(rawcounts) = paste0('iPET_',colnames(rawcounts))
nomcounts =  counts(dds.de, normalized=TRUE)
colnames(nomcounts) = paste0('norm_',colnames(nomcounts))

cormat = cor(nomcounts)
pdf(paste0(cellname, '.', compname, '.heatmapPCC.pdf'), width=6.3, height=6)
pheatmap(cormat, main=paste(cellname, 'PCC DA itx'))
dev.off()

ipetmat  = data.frame(itxID=as.vector(rownames(nomcounts)), rawcounts, nomcounts, stringsAsFactors = FALSE)

#######################################
# Differential Expression Analysis
#######################################
rDE = results(dds.de, name=paste0("condition_", compname), alpha=0.05)
newDE  = lfcShrink(dds.de, coef=paste0("condition_", compname), type="apeglm",  res=rDE)
DEtable = cbind(ipetmat, newDE); #check before merging that rownames are identical

#get annotation
nr = nrow(DEtable)
nodeLR =  strsplit(DEtable$itxID, '_')
nodeL = sapply(nodeLR, '[[',1)
nodeR = sapply(nodeLR, '[[',2)
chL = sapply(strsplit(nodeL, ':'), '[[',1)
chR = sapply(strsplit(nodeR, ':'), '[[',1)
coorL =  sapply(strsplit(nodeL, ':'), '[[',2)
coorR =  sapply(strsplit(nodeR, ':'), '[[',2)
tmpL =strsplit(coorL, '-')
tmpR =strsplit(coorR, '-')
tmpTable = data.frame(chrL = chL, startL = sapply(tmpL, '[[',1), endL=sapply(tmpL, '[[',2) , chrR = chR, startR = sapply(tmpR, '[[',1), endR=sapply(tmpR, '[[',2) , p=0, ipet=1, qval=0) 
tmpname = paste0(cellname, '.tmp.txt')
write.table(tmpTable, sep="\t", quote=F, col.names=F, row.names=F, file=tmpname)
cmd = paste(perlannot, tmpname,  genebed, peakfile, ecDNAfile, roseBed, med1peak,  cellname )
system(cmd)
itx = read.table(paste0(cellname, '.annotated_itx.txt'), header=T, sep="\t", comment.char='', check.names=FALSE, stringsAsFactors=FALSE);
tmpsufs = c('.bedpe.anchors.bed', '.bedpe.anchors.gene.bed', '.cis.sigf.anchors.peak1.intersect.bed', '.cis.sigf.anchors.peak2.intersect.bed','.cis.sigf.anchors.peak3.intersect.bed', '.cis.sigf.anchors.peak4.intersect.bed', '.annotated_itx.txt')
uproL = sapply(strsplit(itx$TSS_L, ";"), function(x) paste(unique(x), collapse = ";")) #unique genes' promoter
uproR = sapply(strsplit(itx$TSS_R, ";"), function(x) paste(unique(x), collapse = ";")) #unique genes' promoter

#SE rank
write.table(cbind(tmpTable[,1:3], 1:nr), sep="\t", quote=F, col.names=F, row.names=F, file=tmpname)
cmd = paste0(bedtools, tmpname, ' -b ', roseBed, ' > ', tmpout)
system(cmd)
o = read.table(tmpout, stringsAsFactors=F)
SErankL = numeric(nr)
SErankL[o$V4] = o$V8

write.table(cbind(tmpTable[,4:6], 1:nr), sep="\t", quote=F, col.names=F, row.names=F, file=tmpname)
cmd = paste0(bedtools, tmpname, ' -b ', roseBed, ' > ', tmpout)
system(cmd)
o = read.table(tmpout, stringsAsFactors=F)
SErankR = numeric(nr)
SErankR[o$V4] = o$V8

mastertab = data.frame(DEtable, ecDNA=itx$ecDNA, PGI=itx$FinalCode, SErank_L=SErankL, SErank_R=SErankR, MED1=itx$MED1, TSS_L=uproL, TSS_R=uproR, chr_L=itx$chrL, start_L=itx$startL, end_L=itx$endL, chr_R=itx$chrR, start_R=itx$startR, end_R=itx$endR)

save.image(paste0(cellname, '.', compname, '.node_itxCis.RData'))
#load(paste0(cellname, '.', compname, '.node_itxCis.RData'))


de.df = as.data.frame(DEtable)
de.df = de.df[order(de.df$padj, abs(de.df$log2FoldChange)),]

thpadj = 0.05
lfc= 1 #log2 scale
signif <- NULL  
bias.cntrl <- NULL  
bias.treat <- NULL  
for (minmean in seq(2,10,by=2)){
    pdedf = subset(de.df, baseMean > minmean & padj >0)
    pdedf$threshold = as.factor(abs(pdedf$log2FoldChange) > lfc & pdedf$padj < thpadj)
    signif = c(signif, table(pdedf$threshold))
    bias.cntrl = c(bias.cntrl, length(which(subset(pdedf, threshold ==T)$log2FoldChange <0)))
    bias.treat= c(bias.treat, length(which(subset(pdedf, threshold ==T)$log2FoldChange >0)))

}

#Just view if iPET affect the up/down regulated numbers
data.frame(NoChange=signif[names(signif)==FALSE], Bias=signif[names(signif)==TRUE], Treat=bias.treat, Control=bias.cntrl)

LRid = strsplit(de.df$itxID, '_')
itxType =  sapply(LRid, function(x){paste0(atype[x[1]],'_', atype[x[2]])})
pairtype = rep('ec-chrom', nrow(de.df))
pairtype[which(itxType == 'ch_ch')] = 'chrom-chrom'
pairtype[which(itxType == 'ec_ec')] = 'ec-ec'
de.df$pair = factor(pairtype)

#####################################################
###  Volcano Plots
#####################################################
chch.col = '#499cbfAA'
ecch.col = '#8e49bf'
ecch.col = 'black'
ecec.col = '#fc03b6c0'

pdedf = de.df; 
colset = c(chch.col, ecch.col, ecec.col)
names(colset) = levels(unique(pdedf$pair))
shpset = c(20,19,23)
names(shpset) = levels(unique(pdedf$pair))
pdedf$color = factor(colset[pdedf$pair])
pdedf$shape = factor(shpset[pdedf$pair])

pdedf$threshold = as.factor(abs(pdedf$log2FoldChange) > lfc & pdedf$padj < thpadj)

a=subset(pdedf, threshold==T)
n.upreg = length(which(a$log2FoldChange > 0)) #46 up drug treatment
n.dnreg = length(which(a$log2FoldChange < 0)) #194 down drug treatment
lpair = table(a$pair)[a$pair]
a$group = factor(paste(names(lpair), lpair))

volc = ggplot(data=pdedf, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(col=threshold), size = 3) +
  scale_color_manual(values=c("#c0c0c0cc","#2942f599")) +
  ggtitle(paste(cellcntrl, n.dnreg, 'vs', celltreat, n.upreg)) + theme_bw(20) + guides(color=guide_legend(title=paste0("log2Fold>",lfc, ";padj<",thpadj)))  +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))

if (length(levels(a$group)) < 3){
    colset = c(chch.col, ecec.col)
    shpset = c(20,18)
}else{
    colset = as.vector(colset)
    shpset = as.vector(shpset)
}
#volca = ggplot(data=a, aes(log2FoldChange, padj)) +
volca = ggplot(data=a, aes(x=log2FoldChange, y=-log10(padj), shape=group)) +
  geom_point(aes(col=group, size=3)) +
  scale_shape_manual(values=shpset) +
  scale_color_manual(values=colset) +
  ggtitle(paste(cellname, 'padj',thpadj, ' FC >', 2^lfc )) + theme_bw(20) + guides(color=guide_legend(title=paste0("Interaction")))  +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))

pdf(paste0(cellname, '.', compname, ".volcanoLoops.pdf"), height =8, width=8)
print(volc) 
print(volca) 
dev.off()

write.table(mastertab, file=paste0(cellname, '.', compname, ".diffDAitxCis.tsv"), sep='\t', quote=F, col.names=T, row.names=F)
system(paste('rm', tmpbed, tmpout, tmpname,  paste0(cellname,tmpsufs, collapse=' ')))
system('rm *.petID.txt')


########################
library(ggplot2)

pdedf$group = factor(ifelse(pdedf$threshold == FALSE, "below threshold 18932",
			ifelse(pdedf$pair == "ec-ec", "ec 59", "non-ec 520")))

pdedf$order = ifelse(pdedf$threshold == FALSE, 3,
			ifelse(pdedf$pair == "ec-ec", 1, 2))


table(pdedf$group)

basic.col = "#c0c0c0cc"
chch.col = "#088dfa99"
ecec.col = "#fd3030c0"
ecch.col = "black"

colset = c(basic.col,ecec.col,chch.col)
names(colset) = levels(unique(pdedf$group))
shpset =c(20,20,20)
names(shpset) = levels(unique(pdedf$group))

    colset = as.vector(colset)
    shpset = as.vector(shpset)

pdedf1 <- pdedf[order(pdedf$order, decreasing=TRUE),]
pdedf = pdedf1
size = 2.5

volc_pdedf = ggplot(data=pdedf, aes(x=log2FoldChange, y=-log10(padj), shape=group)) +
  geom_point(aes(col=group)
                    ,size=size
                    ) +
  #scale_size_manual(values = c(10,10,10))+
  scale_shape_manual(values=shpset) +
  scale_color_manual(values=colset) +
  ggtitle(paste(cellname, 'padj',thpadj, ' FC >', 2^lfc )) + 
  theme_bw(20) + 
  guides(color=guide_legend(title=paste0("Cis interaction")))  +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.border = element_blank(), axis.line = element_line(colour = "black"))

volc_pdedf

ggsave(volc_pdedf, file = paste0(cellname, '.', compname,"_size",size, ".volcanoLoops_merge2.pdf"), height =8, width=8)

########################







