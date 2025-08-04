#given 2 libraries, find differential iPET 
#trans itx: collect all anchors from libraries & replicates
#combine replicas

#merge the overlapping anchors into interaction nodes
#assign id for every node, thus the interaction between 2 nodes will have ID like node1_node2 
#Get a master table for all available samples iPET for each node1_node2 
#normalize with DEseq2
#find overlap with chiadrop anchors

require(data.table)
require(ggplot2)
library(dplyr)
library(pheatmap)
library(ggrepel)

cellname='B171'

cdtsv=paste0('/net/nwgc/vol1/sharing/Wei_Lab/USERS/tungch/tjongh/ec/combine2Reps/',cellname,'.anchors_ecGEMs.annot.tsv')
chiadrop = read.delim(cdtsv)

itxDir = "/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/ChIA-PIPE/diffloop/" #modify accordingly
samplePhenotype = read.table(paste0(itxDir,cellname, "_itxlist.txt")); 
colnames(samplePhenotype) = c("name","condition", "rep") #modify
celltreat <- sub('mock', '',paste(unique(samplePhenotype$condition), collapse='')) #mock must be in samplePhenotype
cellcntrl <- 'mock'
compname <- paste0(celltreat, '_vs_', cellcntrl)
nsample = nrow(samplePhenotype)
comparedLibs  = split(samplePhenotype$name, samplePhenotype$condition)

bedtools='singularity run --bind /net/nwgc/vol1/nobackup/nocleanup/tungch/:/net/nwgc/vol1/nobackup/nocleanup/tungch/ /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools intersect -wo -a '
bedmerge='singularity run --bind /net/nwgc/vol1/nobackup/nocleanup/tungch/:/net/nwgc/vol1/nobackup/nocleanup/tungch/ /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools merge -i '
bedpetools='singularity run --bind /net/nwgc/vol1/nobackup/nocleanup/tungch/:/net/nwgc/vol1/nobackup/nocleanup/tungch/ /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools pairtopair -is -type both -a '

itxlibs = samplePhenotype$name;
itxsuf='.trans.DA.bedpe'

#For annotation

perlannot='/usr/bin/perl /net/nwgc/vol1/sharing/Wei_Lab/USERS/tungch/tjongh/ec/itxAnnotation.4peaks.pl'
genebed='/net/nwgc/vol1/nobackup/nocleanup/tungch/proj-in-situ-chia-pet/annotation/hg38.Genes.itxAnnotation.bed'
regiondir='/net/nwgc/vol1/nobackup/nocleanup/tungch/regions/'

ecDNAfile=paste0(regiondir, cellname, '.ecDNA.bed')
peakfile=paste0(regiondir, cellname, ".Pol2.narrowPeak")
roseBed = paste0(regiondir, cellname , '.SErose.bed')
med1peak = paste0(regiondir, cellname, '.MED1overlap_peak.bed')

message('Annotation files: ', ecDNAfile)
message('Annotation files: ', roseBed )

#-----------------------------------------------------------------------------------------------------------
chnum = 1:23
names(chnum) = paste0('chr',c(1:22,'X'))

#function
getAnchorIndex <- function(itx1, libID){
    #anchorsInfo = read.table(anchorfile, header=F, stringsAsFactors=F)# 6 columns
    #colnames(anchorsInfo) = c('Chr', 'Start', 'End', 'ID', 'num', 'type') 
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
    #out = read.table(tmpout, stringsAsFactors=F)
    numL[out$V4] = out$V9;

    m1 = itx1[,4:6]; m1$idR = 1:n
    write.table(m1, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
    com = paste0(bedtools, tmpbed, ' -b ', anchorfile, ' > ', tmpout)
    system(com)
    out = fread(tmpout, stringsAsFactors=F)
    #out = read.table(tmpout, stringsAsFactors=F)
    numR[out$V4] = out$V9;

    rl = data.frame(L=numL, R=numR)
    rl = subset(rl, L != 0 & R != 0)
    anchorL=apply(rl, 1, min)
    anchorR=apply(rl, 1, max)
    #anchorLab = paste0(anchorsInfo$ID[anchorL], '_', anchorsInfo$ID[anchorR]);
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
itxDir = "/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/ChIA-PIPE/HDtreat/" #modify accordingly
#Reading itx.DA
itxs = list()
for ( i in itxlibs ){
    f1 =  paste0(itxDir, i, "/", i, itxsuf)
    headerset= c('ChrL', 'StartL', 'EndL', 'ChrR', 'StartR', 'EndR', 'iPET')
    message('Read  ', f1)
    #m1 = read.table(f1, header=FALSE, stringsAsFactors=F)
    m1 = fread(f1, header=FALSE, stringsAsFactors=F)
    colnames(m1) = headerset
    itxs[[i]] = m1
}

allitx = bind_rows(itxs) #combine all itx
itx.treat = bind_rows(itxs[comparedLibs[[celltreat]]]) #combine treated
itx.cntrl = bind_rows(itxs[comparedLibs[[cellcntrl]]]) #combine mock

tmpTable=cbind(allitx[,1:6], 0, allitx[,7], 0)
tmpname = paste0(cellname, '.tmp.txt')
write.table(tmpTable, sep="\t", quote=F, col.names=F, row.names=F, file=tmpname)
cmd = paste(perlannot, tmpname,  genebed, peakfile, ecDNAfile, roseBed, med1peak,  cellname )
system(cmd)
itx = read.table(paste0(cellname, '.annotated_itx.txt'), header=T, sep="\t", comment.char='', check.names=FALSE, stringsAsFactors=FALSE);


#selecting ec-chr
#ecchr = subset(itx, ecDNA != '.', select=c( chrL,startL , endL, chrR,   startR,     endR, Nread))

#Let's start with all trans for normalization
ecchr = subset(itx, select=c( chrL,startL , endL, chrR,   startR,     endR, Nread))

#cobine anchor sides
anchors1 = ecchr[,1:3]
anchors2 = ecchr[,4:6]
colnames(anchors1) = colnames(anchors2) = c('CHR','Start','End')
anchors1$num = chnum[anchors1$CHR]
anchors2$num = chnum[anchors2$CHR]


anchors = rbind(anchors1, anchors2)
rm(anchors1,anchors2)

#indexing the combined anchors
tmpbed = paste0(cellname, '.tmp.bed'); 
tmpout = paste0(cellname, '.tmp.out');
write.table(anchors, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0('sort -k4,4n -k2,2n ', tmpbed, ' | ', bedmerge, ' stdin > ', tmpout)  #sort the chnum
system(com)
anchors = fread(tmpout, stringsAsFactors=F) #merging results
#colnames(anchors) =  c('CHR','Start','End'); anchors$ID = 1:nrow(anchors)
colnames(anchors) =  c('CHR','Start','End'); anchors$ID = paste0(anchors$CHR, ':', anchors$Start,'-',anchors$End)
anchors$Num = 1:nrow(anchors)

anchorfile = paste0(cellname, '.anchors.bed')
write.table(anchors, file=anchorfile, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0(bedtools, anchorfile, ' -b ', ecDNAfile, ' > ', tmpout)
system(com)
o = read.table(tmpout, stringsAsFactors=F)
j=which(anchors$ID %in% o$V4) #ecDNA anchors

#annotate anchors for EC or not
atype = rep('ch', nrow(anchors))
atype[j] = 'ec'
names(atype) = anchors$ID
anchors$locus <- atype;
write.table(anchors, file=anchorfile, sep='\t', quote=F, col.names=F, row.names=F) #anchor list with type
anchorsInfo = anchors
colnames(anchorsInfo) = c('Chr', 'Start', 'End', 'ID', 'num', 'type')

#map the itx anchors with ID
anchor.treat = getAnchorIndex(itx.treat, libID=celltreat)
anchor.cntrl = getAnchorIndex(itx.cntrl, libID=cellcntrl)

#alltab = anchor.list %>% reduce(full_join, by="pairID") #merging the tables
alltab = merge(anchor.cntrl, anchor.treat, all=T)
alltab[is.na(alltab)] = 0;
alltab$pairType = paste0(anchorsInfo$type[alltab$numL],'_', anchorsInfo$type[alltab$numR])
iec = which(alltab$pairType != 'ch_ch')
table(alltab$pairType)
#  ch_ch  ch_ec  ec_ch  ec_ec 
# 391832    212  21154      9 

libspec = grep('iPET', colnames(alltab)) #specific column for each lib
grpnames = sub('_iPET', '',colnames(alltab)[libspec])# cntrl and treat names in the order of alltab columns
#Writing input for DEseq2
#for ( i in 1:length(libspec) ){
#    lib1 = grpnames[i]
#    message('Writing labeled itx for ', lib1)
#    write.table(alltab[,c(1,libspec[i])], file=paste0(lib1, '.petID.txt'), sep='\t', quote=F, col.names=F, row.names=F)
#}

countMat = alltab[,libspec] #row is the same order as in alltab, the ec itx are rows iec
colnames(countMat) = grpnames
normFactor = apply(countMat, 2, mean) #median wouldn't work since mock median = 0
normMat = t(apply(countMat, 1, function(x) x/normFactor)) 
colnames(normMat) = paste0('norm_', colnames(normMat))

#normMat0 = normMat
#normMat = countMat
#ipetmat  = data.frame(itxID=as.vector(rownames(nomcounts)), rawcounts, nomcounts, stringsAsFactors = FALSE)
ipetmat = data.frame(alltab[iec,], normMat[iec,], stringsAsFactors = FALSE)
jend = grep('norm', colnames(ipetmat)) #normalized ipet
jR = grep('numR', colnames(ipetmat))
jL = grep('numL', colnames(ipetmat))
ec_ch = subset(ipetmat, pairType == 'ec_ch', select=c(numR, jend))
ch_ec = subset(ipetmat, pairType == 'ch_ec', select=c(numL, jend))

anchors.ec_ch = data.frame(anchorsInfo[ec_ch$numR,1:3], ec_ch[,2:3])
anchors.ch_ec = data.frame(anchorsInfo[ch_ec$numL,1:3], ch_ec[,2:3])
ecConnected.anchors = bind_rows(anchors.ec_ch, anchors.ch_ec)
write.table(ecConnected.anchors, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0('sort -k1,1 -k2,2n ', tmpbed, ' | ', bedmerge, ' stdin -c 4,5 -o sum,sum > ', tmpout)  #sort the chnum
system(com)
# Now tmpout contains chiapet ec-connected chromosomal anchors which will be picked if overlapped with chiadrop


# TSS annotation
ecConnected.nodes = fread(tmpout, stringsAsFactors=F)
colnames(ecConnected.nodes) = colnames(ecConnected.anchors)

tmpTable=cbind(ecConnected.nodes[,1:3],ecConnected.nodes[,1:3], 0, ecConnected.nodes[,4], 0)
tmpname = paste0(cellname, '.tmp.txt')
write.table(tmpTable, sep="\t", quote=F, col.names=F, row.names=F, file=tmpname)
cmd = paste(perlannot, tmpname,  genebed, peakfile, ecDNAfile, roseBed, med1peak,  cellname )
system(cmd)
itx = read.table(paste0(cellname, '.annotated_itx.txt'), header=T, sep="\t", comment.char='', check.names=FALSE, stringsAsFactors=FALSE);
uproL = sapply(strsplit(itx$TSS_L, ";"), function(x) paste(unique(x), collapse = ";")) #unique genes' promoter

ec_chr = data.frame(ecConnected.nodes, ecDNA=itx$ecDNA, TSS_L=uproL)



#Overlap chiapet with chiadrop
chiadrop.P = subset(chiadrop
                    , ec == FALSE & PGI =='P'
                    , select=c(Chr, Start, End, nGEMs,MergeFragID, TSS)
                    )

colnames(chiadrop)

#counts
N.chiadrop <- nrow(chiadrop)
message("N.chiadrop: ", N.chiadrop)
# N.chiadrop: 3838
N.chiadrop.P <- nrow(chiadrop.P)
message("N.chiadrop.P: ", N.chiadrop.P)
#N.chiadrop.P: 2751

tmpout <- fread(tmpout, stringsAsFactors=F)
N.tmpout <- nrow(tmpout)
message("N.tmpout: ", N.tmpout)
#N.tmpout: 11980

tmpbed = paste0(cellname, '.tmp.bed'); 
tmpout = paste0(cellname, '.tmp.out');
colnames(chiadrop.P) = paste0('chiadrop_',colnames(chiadrop.P))
write.table(chiadrop.P, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)



com = paste0(bedtools, tmpout, ' -b ', tmpbed, ' > ', tmpname)
system(com)
oo  = fread(tmpname, stringsAsFactors=F) #results of chiadrop overlap

#counts
N.oo <- nrow(oo)
message("N.oo: ", N.oo)
#N.oo: 1662

colnames(oo) = c(colnames(ecConnected.anchors), colnames(chiadrop.P), 'bpOverlap')



oo$chr_anchor <- paste0(oo$Chr,":",oo$Start,"-",oo$End)

head(oo)

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

sum(merged_ecConnected$norm_mock != 0)
sum(merged_ecConnected$norm_drug != 0)
# > sum(merged_ecConnected$norm_mock != 0)
# [1] 945
# > sum(merged_ecConnected$norm_drug != 0)
# [1] 1234


## unique gene number check
# separate
library(tidyr)
colnames(merged_ecConnected)
outtab <-  separate_rows(merged_ecConnected, chiadrop_TSS, sep = ";")
outtab = as.data.frame(outtab)
dim(outtab)
# [1] 2362   13
table(duplicated(outtab$chiadrop_TSS))
# FALSE  TRUE 
#  2198   164 
length(unique(outtab$chiadrop_TSS))
# [1] 2198

#

Values_norm_mock = as.numeric(merged_ecConnected$norm_mock[merged_ecConnected$norm_mock != 0])
Values_norm_drug = as.numeric(merged_ecConnected$norm_drug[merged_ecConnected$norm_drug != 0])

wilcox.test(Values_norm_drug, 
            Values_norm_mock, 
            alternative='less')$p.value
# 9.474219e-47

nA_v = length(Values_norm_mock)
# [1] 945


nA_1 = length(Values_norm_drug)
# [1] 1234


cplabel = c(rep("mock", nA_v), rep("HD", nA_1))
cp.df = data.frame(normFreq=c(Values_norm_mock, Values_norm_drug),  cond=factor(cplabel, levels=c("mock","HD"))) #check if grpnames call is consistent

summary(Values_norm_mock)
summary(Values_norm_drug)

# > summary(Values_norm_mock)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.356   2.356   2.356   4.341   4.712  27.095 
# > summary(Values_norm_drug)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.542   1.542   2.313   3.328   4.626  20.045

##Interaction violin
library(ggplot2)
library(ggbreak) 
library(scales)
library(gg.gap)
library(ggbeeswarm)

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

  # geom_beeswarm() +
  # geom_quasirandom(dodge.width=.75, col=2) +
  # geom_point(position=position_jitterdodge(jitter.width = 0,dodge.width = 2)) +
  #scale_color_manual(values=c("#B099BF", "#8FBF88")) +
  scale_fill_manual(values=c("mock" = "#7fbf7b", "HD" = "#af8dc3")) +
  geom_boxplot(width=0.03, color="black", alpha=10) +
  # geom_boxplot(width=0.15, color="black", alpha=0.2, outlier.shape = NA) +
  # ylim(0, 50) +
  expand_limits(y=c(1, 100)) +
  scale_y_continuous(
    # labels = scientific,
    trans='log10',
    # trans = pseudo_log_trans(base = 10),
    breaks = c(1, 10, 100),
    # minor_breaks = mb,
    guide = "axis_logticks",
    # labels = label_number(accuracy = 1),
    # minor_breaks = log10_minor_breaks()
    ) +
        #  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
        #       labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
          plot_theme
dev.off()

# pdf(paste0(cellname, '.', compname, '.diff_ecTrans_violin_width6.pdf'), width=6, height=6)
# ggplot(a, aes(x=cond, y=normFreq, fill=cond)) + 
#   geom_violin(width=0.7, bw =0.05) +
#   scale_x_discrete(labels=c(paste0(unique(cplabel)[1]," (n = ",nA_v,")"), paste0(unique(cplabel)[2]," (n = ",nA_1,")"))) +
#   xlab("") +
#   #scale_color_manual(values=c("#B099BF", "#8FBF88")) +
#   scale_fill_manual(values=c("mock" = "#7fbf7b", "HD" = "#af8dc3")) +
#   geom_boxplot(width=0.03, color="black", alpha=10) +
#   # geom_boxplot(width=0.15, color="black", alpha=0.2, outlier.shape = NA) +
#   #ylim(0, 120) +
#   scale_y_continuous(
#     trans = pseudo_log_trans(base = 10),
#     breaks = c(1, 10, 100),
#     # minor_breaks = mb,
#     guide = "axis_logticks",
#     labels = label_number(accuracy = 1),
#     # minor_breaks = log10_minor_breaks()
#     ) +
#   theme(legend.position="bottom",
#         panel.spacing.x = unit(1.5, "lines")) +
#           plot_theme
# dev.off()



# write master table
write.table(merged_ecConnected  , file=paste0(cellname, '.', compname, '.CombinedReplicate.Trans.chiadrop_genes_nonPairwise.tsv'), sep='\t', quote=F, col.names=T, row.names=F)


save.image(paste0(cellname, '.', compname, '.CombinedReplicate.Trans.chiadrop_genes_nonPairwise.RData'))


# stats
  sampleList = c(paste0("mock-1"),
                paste0("mock-2"),
                paste0("HD-1"),
                paste0("HD-2"))
  alltab_id = c(4,4,5,5)

  results <- data.frame(LIBname = character(),
                        Sample = character(),
                        total_trans.DA = numeric(),
                        total_node.based_trans_chrMclean = numeric(),

                        ecTrans = numeric(),
                        ecTrans_TSS = numeric(),
                        ecTrans_chiadrop_TSS = numeric(),                        

                        ecTrans_chiadrop_uniqueGenes = numeric()
                        )


  for (j in 1:4) {

    results[nrow(results) + 1, ] <- c(names(itxs[j]),
                                        sampleList[j],
                                        nrow(as.data.frame(itxs[j])),
                                        sum(as.data.frame(alltab)[,alltab_id[j]] != 0),

                                        sum(as.data.frame(ecConnected.anchors)[,alltab_id[j]] != 0),
                                        sum(ec_chr$TSS_L != "." & ec_chr[,alltab_id[j]] != 0),
                                        length(merged_ecConnected$chiadrop_MergeFragID[as.data.frame(merged_ecConnected)[,alltab_id[j]+1] != 0]),                                        

                                        length(unique(outtab$chiadrop_TSS[as.data.frame(outtab)[,alltab_id[j]+1] != 0]))
                                    )
    

  }

    system("rm *tmp* *_itx.txt *petID* B171*bed")

out = t(results)
write.table(out, file=paste0(cellname,'.',compname,'.ecTrans_chiadrop_based.stats.txt'), quote=F, row.names=T, col.names=F, sep="\t")
