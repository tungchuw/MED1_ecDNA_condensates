require(data.table)
require(ggplot2)
library(dplyr)
library(gg.gap)
library(ggbeeswarm)

cellname='PC3DM+'

cdtsv=paste0('./',cellname,'.anchors_ecGEMs.annot.tsv')
chiadrop = read.delim(cdtsv)

itxDir = "./ChIA-PIPE/diffloop/" #modify accordingly
samplePhenotype = read.table(paste0(itxDir,cellname, "_itxlist.txt"));
colnames(samplePhenotype) = c("name","condition", "rep") #modify
celltreat <- sub('mock', '',paste(unique(samplePhenotype$condition), collapse='')) #mock must be in samplePhenotype
cellcntrl <- 'mock'
compname <- paste0(celltreat, '_vs_', cellcntrl)
nsample = nrow(samplePhenotype)
comparedLibs  = split(samplePhenotype$name, samplePhenotype$condition)

bedtools='singularity run ./singularity/chipseqtools.sif bedtools intersect -wo -a '
bedmerge='singularity run ./singularity/chipseqtools.sif bedtools merge -i '
bedpetools='singularity run ./singularity/chipseqtools.sif bedtools pairtopair -is -type both -a '

itxlibs = samplePhenotype$name;
itxsuf='.trans.DA.bedpe'

#For annotation
perlannot='/usr/bin/perl ./itxAnnotation.4peaks.pl'
genebed='./hg38.Genes.itxAnnotation.bed'
regiondir='./'

ecDNAfile=paste0(regiondir, cellname, '.ecDNA.bed')
peakfile=paste0(regiondir, cellname, ".Pol2.narrowPeak")
roseBed = paste0(regiondir, cellname , '.SErose.bed')
med1peak = paste0(regiondir, cellname, '.MED1overlap_peak.bed')

#-----------------------------------------------------------------------------------------------------------
chnum = 1:23
names(chnum) = paste0('chr',c(1:22,'X'))

#function
getAnchorIndex <- function(itx1, libID){
    n = nrow(itx1)
    numL = numeric(n);
    numR = numeric(n);
    names(numL) = names(numR) = 1:n
    tmpbed = 'tmp.bed'; tmpout = 'tmp.txt'

    m1 = itx1[,1:3]; m1$idL = 1:n
    write.table(m1, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
    com = paste0(bedtools, tmpbed, ' -b ', anchorfile, ' > ', tmpout)
    system(com)
    out = fread(tmpout, stringsAsFactors=F)
    numL[out$V4] = out$V9;

    m1 = itx1[,4:6]; m1$idR = 1:n
    write.table(m1, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
    com = paste0(bedtools, tmpbed, ' -b ', anchorfile, ' > ', tmpout)
    system(com)
    out = fread(tmpout, stringsAsFactors=F)
    numR[out$V4] = out$V9;

    rl = data.frame(L=numL, R=numR)
    rl = subset(rl, L != 0 & R != 0)
    anchorL=apply(rl, 1, min)
    anchorR=apply(rl, 1, max)
    anchorLab = paste0(anchorL, '_', anchorR);
    idPETs = data.frame(pairID=anchorLab, iPET=itx1[,7], stringsAsFactors=FALSE)
    idPETs = aggregate(iPET~pairID, idPETs, sum)
    colnames(idPETs)[2] = paste0(libID, '_', colnames(idPETs)[2])
    #recover the region by numeric id
    s = strsplit(idPETs$pairID,'_')
    numL =  as.numeric(sapply(s, '[[',1 ));
    numR =  as.numeric(sapply(s, '[[',2 ));
    idPETs$numL = numL;
    idPETs$numR = numR;
    system('rm tmp.bed tmp.txt')
    return(idPETs)
}

#==============================
itxDir = "./" #modify accordingly
itxs = list()
for ( i in itxlibs ){
    f1 =  paste0(itxDir, i, "/", i, itxsuf)
    headerset= c('ChrL', 'StartL', 'EndL', 'ChrR', 'StartR', 'EndR', 'iPET')
    m1 = fread(f1, header=FALSE, stringsAsFactors=F)
    colnames(m1) = headerset
    itxs[[i]] = m1
}

allitx = bind_rows(itxs)
itx.treat = bind_rows(itxs[comparedLibs[[celltreat]]])
itx.cntrl = bind_rows(itxs[comparedLibs[[cellcntrl]]])

tmpTable=cbind(allitx[,1:6], 0, allitx[,7], 0)
tmpname = paste0(cellname, '.tmp.txt')
write.table(tmpTable, sep="\t", quote=F, col.names=F, row.names=F, file=tmpname)
cmd = paste(perlannot, tmpname,  genebed, peakfile, ecDNAfile, roseBed, med1peak,  cellname )
system(cmd)
itx = read.table(paste0(cellname, '.annotated_itx.txt'), header=T, sep="\t", comment.char='', check.names=FALSE, stringsAsFactors=FALSE);

ecchr = subset(itx, select=c( chrL,startL , endL, chrR,   startR,     endR, Nread))

anchors1 = ecchr[,1:3]
anchors2 = ecchr[,4:6]
colnames(anchors1) = colnames(anchors2) = c('CHR','Start','End')
anchors1$num = chnum[anchors1$CHR]
anchors2$num = chnum[anchors2$CHR]

anchors = rbind(anchors1, anchors2)
rm(anchors1,anchors2)

tmpbed = paste0(cellname, '.tmp.bed');
tmpout = paste0(cellname, '.tmp.out');
write.table(anchors, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0('sort -k4,4n -k2,2n ', tmpbed, ' | ', bedmerge, ' stdin > ', tmpout)
system(com)
anchors = fread(tmpout, stringsAsFactors=F)
colnames(anchors) =  c('CHR','Start','End'); anchors$ID = paste0(anchors$CHR, ':', anchors$Start,'-',anchors$End)
anchors$Num = 1:nrow(anchors)

anchorfile = paste0(cellname, '.anchors.bed')
write.table(anchors, file=anchorfile, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0(bedtools, anchorfile, ' -b ', ecDNAfile, ' > ', tmpout)
system(com)
o = read.table(tmpout, stringsAsFactors=F)
j=which(anchors$ID %in% o$V4)

atype = rep('ch', nrow(anchors))
atype[j] = 'ec'
names(atype) = anchors$ID
anchors$locus <- atype;
write.table(anchors, file=anchorfile, sep='\t', quote=F, col.names=F, row.names=F)
anchorsInfo = anchors
colnames(anchorsInfo) = c('Chr', 'Start', 'End', 'ID', 'num', 'type')

anchor.treat = getAnchorIndex(itx.treat, libID=celltreat)
anchor.cntrl = getAnchorIndex(itx.cntrl, libID=cellcntrl)

alltab = merge(anchor.cntrl, anchor.treat, all=T)
alltab[is.na(alltab)] = 0;
alltab$pairType = paste0(anchorsInfo$type[alltab$numL],'_', anchorsInfo$type[alltab$numR])
iec = which(alltab$pairType != 'ch_ch')

libspec = grep('iPET', colnames(alltab))
grpnames = sub('_iPET', '',colnames(alltab)[libspec])

countMat = alltab[,libspec]
colnames(countMat) = grpnames
normFactor = apply(countMat, 2, mean)
normMat = t(apply(countMat, 1, function(x) x/normFactor))
colnames(normMat) = paste0('norm_', colnames(normMat))

ipetmat = data.frame(alltab[iec,], normMat[iec,], stringsAsFactors = FALSE)
jend = grep('norm', colnames(ipetmat))
jR = grep('numR', colnames(ipetmat))
jL = grep('numL', colnames(ipetmat))
ec_ch = subset(ipetmat, pairType == 'ec_ch', select=c(numR, jend))
ch_ec = subset(ipetmat, pairType == 'ch_ec', select=c(numL, jend))

anchors.ec_ch = data.frame(anchorsInfo[ec_ch$numR,1:3], ec_ch[,2:3])
anchors.ch_ec = data.frame(anchorsInfo[ch_ec$numL,1:3], ch_ec[,2:3])
ecConnected.anchors = bind_rows(anchors.ec_ch, anchors.ch_ec)
write.table(ecConnected.anchors, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0('sort -k1,1 -k2,2n ', tmpbed, ' | ', bedmerge, ' stdin -c 4,5 -o sum,sum > ', tmpout)
system(com)

ecConnected.nodes = fread(tmpout, stringsAsFactors=F)
colnames(ecConnected.nodes) = colnames(ecConnected.anchors)

tmpTable=cbind(ecConnected.nodes[,1:3],ecConnected.nodes[,1:3], 0, ecConnected.nodes[,4], 0)
tmpname = paste0(cellname, '.tmp.txt')
write.table(tmpTable, sep="\t", quote=F, col.names=F, row.names=F, file=tmpname)
cmd = paste(perlannot, tmpname,  genebed, peakfile, ecDNAfile, roseBed, med1peak,  cellname )
system(cmd)
itx = read.table(paste0(cellname, '.annotated_itx.txt'), header=T, sep="\t", comment.char='', check.names=FALSE, stringsAsFactors=FALSE);
uproL = sapply(strsplit(itx$TSS_L, ";"), function(x) paste(unique(x), collapse = ";"))

ec_chr = data.frame(ecConnected.nodes, ecDNA=itx$ecDNA, TSS_L=uproL)

chiadrop.P = subset(chiadrop
                    , ec == FALSE & PGI =='P'
                    , select=c(Chr, Start, End, nGEMs,MergeFragID, TSS)
                    )

tmpout_df <- fread(tmpout, stringsAsFactors=F) # Re-reading the merged anchors file

tmpbed = paste0(cellname, '.tmp.bed');
tmpout_overlap = paste0(cellname, '.tmp.out'); # Use a different name to avoid confusion
colnames(chiadrop.P) = paste0('chiadrop_',colnames(chiadrop.P))
write.table(chiadrop.P, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
write.table(tmpout_df, file=tmpout, sep='\t', quote=F, col.names=F, row.names=F) # Write the correct data for intersection

com = paste0(bedtools, tmpout, ' -b ', tmpbed, ' > ', tmpname)
system(com)
oo  = fread(tmpname, stringsAsFactors=F)

colnames(oo) = c(colnames(ecConnected.anchors), colnames(chiadrop.P), 'bpOverlap')
oo$chr_anchor <- paste0(oo$Chr,":",oo$Start,"-",oo$End)

merged_ecConnected <- oo %>%
                      group_by(chr_anchor) %>%
                      summarise(
                        Chr = paste(unique(Chr), collapse = ";"),
                        Start = paste(unique(Start), collapse = ";"),
                        End = paste(unique(End), collapse = ";"),
                        norm_mock = paste(unique(norm_mock), collapse = ";"),
                        norm_drug = paste(unique(norm_drug), collapse = ";"),
                        chiadrop_Chr = paste(unique(chiadrop_Chr), collapse = ";"),
                        chiadrop_Start = paste(unique(chiadrop_Start), collapse = ";"),
                        chiadrop_End = paste(unique(chiadrop_End), collapse = ";"),
                        chiadrop_nGEMs = paste(chiadrop_nGEMs, collapse = ";"),
                        chiadrop_MergeFragID = paste(unique(chiadrop_MergeFragID), collapse = ";"),
                        chiadrop_TSS = paste(unique(chiadrop_TSS), collapse = ";"),
                        bpOverlap = paste(unique(bpOverlap), collapse = ";"),
                        .groups = 'drop'
                      )

merged_ecConnected = as.data.frame(merged_ecConnected)

Values_norm_mock = as.numeric(merged_ecConnected$norm_mock[merged_ecConnected$norm_mock != 0])
Values_norm_drug = as.numeric(merged_ecConnected$norm_drug[merged_ecConnected$norm_drug != 0])

wilcox.test(Values_norm_drug,
            Values_norm_mock,
            alternative='less')

nA_v = length(Values_norm_mock)
nA_1 = length(Values_norm_drug)

cplabel = c(rep("mock", nA_v), rep("HD", nA_1))
cp.df = data.frame(normFreq=c(Values_norm_mock, Values_norm_drug),  cond=factor(cplabel, levels=c("mock","HD")))

##Interaction violin
library(ggplot2)
library(scales)

a=subset(cp.df)

plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 12),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.border = element_blank(),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

pdf(paste0(cellname, '.', compname, '.CombinedReplicate.Trans.chiadrop_genes_nonPairwise_min1-max100.pdf'), width=6, height=6)
ggplot(a, aes(x=cond, y=normFreq, fill=cond)) +
  geom_violin(width=0.7, bw =0.07) +
  scale_x_discrete(labels=c(paste0(unique(cplabel)[1]," (n = ",nA_v,")"), paste0(unique(cplabel)[2]," (n = ",nA_1,")"))) +
  xlab("") +
  scale_fill_manual(values=c("mock" = "#7fbf7b", "HD" = "#af8dc3")) +
  geom_boxplot(width=0.03, color="black", alpha=10) +
  expand_limits(y=c(1, 100)) +
  scale_y_continuous(
    labels = scientific,
    trans='log10',
    breaks = c(0, 1, 10, 100),
    guide = "axis_logticks",
    ) +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
          plot_theme
dev.off()

# write master table
write.table(merged_ecConnected  , file=paste0(cellname, '.', compname, '.CombinedReplicate.Trans.chiadrop_genes_nonPairwise.tsv'), sep='\t', quote=F, col.names=T, row.names=F)

system("rm *tmp* *_itx.txt *petID* PC3DM+*bed")
