# Fig.5B_XX.R - Script to analyze ecDNA-connected chromosomal interactions and overlaps with ChIA-Drop

# Load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(DESeq2)
  library(dplyr)
  library(pheatmap)
  library(ggrepel)
  library(scales)
  library(tidyverse)
  library(tidyr)
  library(ggbreak)
  library(gg.gap)
  library(ggbeeswarm)
})

# === CONFIGURATION ===
cellname <- "PC3DM+"
itxsuf <- '.e500.clusters.trans.chiasig.gz'

# Directory paths
base_dir <- "./data/misc/"
cell_dir <- "./data/misc/PC3DM/"
script_dir <- "./scripts/"

# Input files
cdtsv <- paste0('./data/misc/PC3DM/chiadrop/', cellname, '.anchors_ecGEMs.annot.tsv')
chiadrop = read.delim(cdtsv)
sample_file <- file.path(cell_dir, paste0("chiatac/",cellname, "_all_itxlist.txt"))
samplePhenotype <- read.table(sample_file)
colnames(samplePhenotype) <- c("name","condition", "rep")
itxlibs <- samplePhenotype$name

# Bedtools command wrappers
base_sif <- './singularity/chipseqtools.sif '
bedcmd <- function(tool) paste("singularity run ", base_sif, tool)
bedtools <- bedcmd("bedtools intersect -wo -a ")
bedmerge <- bedcmd("bedtools merge -i ")
bedpetools <- bedcmd("bedtools pairtopair -is -type both -a ")

# Annotation files
ecDNAfile <- file.path(cell_dir, paste0(cellname, ".ecDNA.bed"))
peakfile <- file.path(cell_dir, paste0(cellname, ".Pol2.narrowPeak"))
roseBed <- file.path(cell_dir, paste0(cellname , ".SErose.bed"))
med1peak <- file.path(cell_dir, paste0(cellname, ".MED1overlap_peak.bed"))
genebed <- file.path(base_dir, "hg38.Genes.itxAnnotation.bed")
perlannot <- paste0("/usr/bin/perl ", file.path(script_dir, "itxAnnotation.4peaks.pl"))

# Output file paths
tmpbed <- paste0(cellname, '.tmp.bed')
tmpout <- paste0(cellname, '.tmp.out')
tmpname <- paste0(cellname, '.tmp.txt')
anchorfile <- paste0(cellname, '.anchors.bed')


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
    anchorL[out$V4] = out$V8;

    m1 = itx1[,4:6]; m1$idR = 1:n
    write.table(m1, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
    com = paste0(bedtools, tmpbed, ' -b ', anchorfile, ' > ', tmpout)
    system(com)
    out = fread(tmpout, stringsAsFactors=F)
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
# Reading data
itxs = list()
for ( i in itxlibs ){
    f1 =  paste0(cell_dir, "/chiatac/", i, itxsuf)
    headerset= c('ChrL', 'StartL', 'EndL', 'ChrR', 'StartR', 'EndR', 'iPET')
    message('Read  ', f1)
    m1 = fread(f1, header=FALSE, stringsAsFactors=F)[,1:7]
    colnames(m1) = headerset
    itxs[[i]] = m1
}

# Get anchors
allitx = bind_rows(itxs) #combine all itx
anchors1 = allitx[,1:3]
anchors2 = allitx[,4:6]
colnames(anchors1) = colnames(anchors2) = c('CHR','Start','End')
anchors = rbind(anchors1, anchors2)
rm(anchors1,anchors2)

# Indexing the combined anchors
tmpbed = paste0(cellname, '.tmp.bed'); 
tmpout = paste0(cellname, '.tmp.out');
write.table(anchors, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0('sort -k1,1 -k2,2n ', tmpbed, ' | ', bedmerge, ' stdin > ', tmpout) 
system(com)
anchors = fread(tmpout, stringsAsFactors=F)
colnames(anchors) =  c('CHR','Start','End'); anchors$ID = paste0(anchors$CHR, ':', anchors$Start,'-',anchors$End)
anchorfile = paste0(cellname, '.anchors.bed')
write.table(anchors, file=anchorfile, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0(bedtools, anchorfile, ' -b ', ecDNAfile, ' > ', tmpout)
system(com)
o = read.table(tmpout, stringsAsFactors=F)
atype = rep('ch', nrow(anchors))
atype[which(anchors$ID %in% o$V4)] = 'ec'
names(atype) = anchors$ID
anchors$locus <- atype;

# Map the itx anchors with ID
anchor.list = list()
for ( i in itxlibs ){
    message('Getting pairID for ', i)
    anchor.list[[i]] = getAnchorIndex(itxs[[i]], libID=i)
}

# Merge to a matrix
alltab = merge(anchor.list[[1]],anchor.list[[2]], all=T)
alltab_2 = merge(alltab,anchor.list[[3]], all=T)

# Separate node_L and node_R
alltab_3 <- alltab_2 %>% 
  mutate(pair_L = sub("_.*", "", pairID),
         pair_R = sub(".*_", "", pairID)) %>%
  separate(pair_L, into = c("chr_L", "range_L"), sep = ":") %>% 
  separate(range_L, into = c("start_L", "end_L"), sep = "-") %>%
  separate(pair_R, into = c("chr_R", "range_R"), sep = ":") %>% 
  separate(range_R, into = c("start_R", "end_R"), sep = "-")

# Filter chrM
rawcounts = subset(alltab_3, chr_L != "chrM" & chr_R != "chrM", select = c(1:4))

# Normalization
rawcounts[is.na(rawcounts)] = 0;
row.names(rawcounts) = rawcounts$pairID
rawcounts = rawcounts[,2:4]
normFactor = apply(rawcounts, 2, mean) 
normMat = t(apply(rawcounts, 1, function(x) x/normFactor)) 
colnames(normMat) = paste0('norm_', colnames(normMat))
DEtable = cbind(rawcounts, normMat); 
DEtable$itxID = row.names(DEtable)
head(DEtable)

# Get annotation
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

mastertab_filtered = mastertab
ec_ch = subset(mastertab_filtered, ecDNA == "L", select = c(7,1:6,10,13:20))
ch_ec = subset(mastertab_filtered, ecDNA == "R", select = c(7,1:6,11,13:20))
ch_ec$new_pairID = paste0(ch_ec$chr_R, ":" ,ch_ec$start_R,'-', ch_ec$end_R,"_",ch_ec$chr_L, ":" ,ch_ec$start_L,'-', ch_ec$end_L)

new.ch_ec = data.frame(itxID = ch_ec$new_pairID,
                    AT4.40.19m20.vector_iPET = ch_ec$AT4.40.19m20.vector_iPET,
                    AT4.40.1m2.SE1_iPET = ch_ec$AT4.40.1m2.SE1_iPET,
                    AT4.40.9m10.combo_SE_iPET = ch_ec$AT4.40.9m10.combo_SE_iPET,
                    norm_AT4.40.19m20.vector_iPET = ch_ec$norm_AT4.40.19m20.vector_iPET,
                    norm_AT4.40.1m2.SE1_iPET = ch_ec$norm_AT4.40.1m2.SE1_iPET,
                    norm_AT4.40.9m10.combo_SE_iPET = ch_ec$norm_AT4.40.9m10.combo_SE_iPET,
                    SErank_L = ch_ec$SErank_R,
                      TSS_L = ch_ec$TSS_R,
                      TSS_R = ch_ec$TSS_L,
                      chr_L = ch_ec$chr_R,
                      start_L = ch_ec$start_R,
                      end_L = ch_ec$end_R,
                      chr_R = ch_ec$chr_L,
                      start_R = ch_ec$start_L,
                      end_R = ch_ec$end_L
)
ecConnected.info = bind_rows(ec_ch, new.ch_ec)

# Merge the same anchors and add up the norm_counts
merged_ecConnected <- ecConnected.info %>%
                      group_by(itxID) %>%
                      summarise(
                        AT4.40.19m20.vector_iPET = paste(unique(AT4.40.19m20.vector_iPET), collapse = ";"),
                        AT4.40.1m2.SE1_iPET = paste(unique(AT4.40.1m2.SE1_iPET), collapse = ";"),
                        AT4.40.9m10.combo_SE_iPET = paste(unique(AT4.40.9m10.combo_SE_iPET), collapse = ";"),
                        norm_AT4.40.19m20.vector_iPET = paste(unique(norm_AT4.40.19m20.vector_iPET), collapse = ";"),
                        norm_AT4.40.1m2.SE1_iPET = paste(unique(norm_AT4.40.1m2.SE1_iPET), collapse = ";"),
                        norm_AT4.40.9m10.combo_SE_iPET = paste(unique(norm_AT4.40.9m10.combo_SE_iPET), collapse = ";"),
                        TSS_L = paste(unique(TSS_L), collapse = ";"),
                        TSS_R = paste(unique(TSS_R), collapse = ";"),
                        SErank_L = paste(unique(SErank_L), collapse = ";"),
                        chr_L = paste(unique(chr_L), collapse = ";"),
                        start_L = paste(unique(start_L), collapse = ";"),
                        end_L = paste(unique(end_L), collapse = ";"),
                        chr_R = paste(unique(chr_R), collapse = ";"),
                        start_R = paste(unique(start_R), collapse = ";"),
                        end_R = paste(unique(end_R), collapse = ";"),                        
                        .groups = 'drop' 
                      )

merged_ecConnected = as.data.frame(merged_ecConnected)

# Made table for chiadrop
ecConnected_regions <- merged_ecConnected[,c(14:16,11:13,1:10)]
write.table(ecConnected_regions, file=tmpout, sep='\t', quote=F, col.names=F, row.names=F)
# Now tmpout contains chiapet ec-connected chromosomal anchors which will be picked if overlapped with chiadrop

# Overlap chiapet with chiadrop
chiadrop.P = subset(chiadrop
                    , ec == FALSE & PGI =='P'
                    , select=c(Chr, Start, End, nGEMs,MergeFragID, TSS)
                    )

#counts
N.chiadrop <- nrow(chiadrop)
message("N.chiadrop: ", N.chiadrop)
N.chiadrop.P <- nrow(chiadrop.P)
message("N.chiadrop.P: ", N.chiadrop.P)
tmpout <- fread(tmpout, stringsAsFactors=F)
N.tmpout <- nrow(tmpout)
message("N.tmpout: ", N.tmpout)
tmpbed = paste0(cellname, '.tmp.bed'); 
tmpout = paste0(cellname, '.tmp.out');
colnames(chiadrop.P) = paste0('chiadrop_',colnames(chiadrop.P))
write.table(chiadrop.P, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0(bedtools, tmpout, ' -b ', tmpbed, ' > ', tmpname)
system(com)
oo  = fread(tmpname, stringsAsFactors=F) #results of chiadrop overlap
N.oo <- nrow(oo)
message("N.oo: ", N.oo)
colnames(oo) = c(colnames(ecConnected_regions), colnames(chiadrop.P), 'bpOverlap')

# Reomve duplicated itxID
oo_unique <- oo %>%
            group_by(itxID) %>%
            summarise(
              Chr_R = paste(unique(chr_R), collapse = ";"),
              Start_R = paste(unique(start_R), collapse = ";"),
              End_R = paste(unique(end_R), collapse = ";"),
              Chr_L = paste(unique(chr_L), collapse = ";"),
              Start_L = paste(unique(start_L), collapse = ";"),
              End_L = paste(unique(end_L), collapse = ";"),
              AT4.40.19m20.vector_iPET = paste(unique(AT4.40.19m20.vector_iPET), collapse = ";"),
              AT4.40.1m2.SE1_iPET = paste(unique(AT4.40.1m2.SE1_iPET), collapse = ";"),
              AT4.40.9m10.combo_SE_iPET = paste(unique(AT4.40.9m10.combo_SE_iPET), collapse = ";"),
              norm_AT4.40.19m20.vector_iPET = paste(unique(norm_AT4.40.19m20.vector_iPET), collapse = ";"),
              norm_AT4.40.1m2.SE1_iPET = paste(unique(norm_AT4.40.1m2.SE1_iPET), collapse = ";"),
              norm_AT4.40.9m10.combo_SE_iPET = paste(unique(norm_AT4.40.9m10.combo_SE_iPET), collapse = ";"),
              TSS_L = paste(unique(TSS_L), collapse = ";"),
              TSS_R = paste(unique(TSS_R), collapse = ";"),
              SErank_L = paste(unique(SErank_L), collapse = ";"),

              chiadrop_Chr = paste(chiadrop_Chr, collapse = ";"),
              chiadrop_Start = paste(chiadrop_Start, collapse = ";"),
              chiadrop_End = paste(chiadrop_End, collapse = ";"),
              chiadrop_nGEMs = paste(chiadrop_nGEMs, collapse = ";"),
              chiadrop_MergeFragID = paste(chiadrop_MergeFragID, collapse = ";"),
              chiadrop_TSS = paste(chiadrop_TSS, collapse = ";"),
              bpOverlap = paste(bpOverlap, collapse = ";"),
              .groups = 'drop' 
            )

outtab_id3_TSS <- as.data.frame(oo_unique)
outtab_id3_TSS_SE1 = subset(outtab_id3_TSS, SErank_L == "1")

Values_morm_vector_19n20 = as.numeric(outtab_id3_TSS_SE1$norm_AT4.40.19m20.vector_iPET[outtab_id3_TSS_SE1$norm_AT4.40.19m20.vector_iPET != 0])
Values_morm_SE1_1n2 = as.numeric(outtab_id3_TSS_SE1$norm_AT4.40.1m2.SE1_iPET[outtab_id3_TSS_SE1$norm_AT4.40.1m2.SE1_iPET != 0])
Values_morm_comboSE_9n10 = as.numeric(outtab_id3_TSS_SE1$norm_AT4.40.9m10.combo_SE_iPET[outtab_id3_TSS_SE1$norm_AT4.40.9m10.combo_SE_iPET != 0])

wilcox.test(Values_morm_SE1_1n2, 
            Values_morm_vector_19n20, 
            alternative='less')$p.value
# p-value = 1.369851e-158

wilcox.test(Values_morm_comboSE_9n10, 
            Values_morm_vector_19n20, 
            alternative='less')$p.value
# p-value = 1.640885e-175

nA_v = length(Values_morm_vector_19n20)
# [1] 578
nA_1 = length(Values_morm_SE1_1n2)
# [1] 783
nA_c = length(Values_morm_comboSE_9n10)
# [1] 904

cplabel = c(rep("vector", nA_v), rep("SE1", nA_1), rep("combo_SE", nA_c))
cp.df = data.frame(normFreq=c(Values_morm_vector_19n20, Values_morm_SE1_1n2, Values_morm_comboSE_9n10),  cond=factor(cplabel, levels=c("vector","SE1","combo_SE"))) #check if grpnames call is consistent

# Violin
a=subset(cp.df)
pdf(paste0(cellname, '.MergedFastq.Trans.chiadrop_genes_nonPairwise_SE1filtered_bedtools.pdf'), width=6, height=6)
ggplot(a, aes(x=cond, y=normFreq, fill=cond)) + 
  geom_violin(width=0.7, bw =0.05) +
  scale_x_discrete(labels=c(paste0(unique(cplabel)[1]," (n = ",nA_v,")"), paste0(unique(cplabel)[2]," (n = ",nA_1,")"), paste0(unique(cplabel)[3]," (n = ",nA_c,")"))) +
  xlab("") +
  scale_fill_manual(values=c("vector" = "#7fbf7b", "SE1" = "#af8dc3", "combo_SE" = "#af8dc3")) +
  geom_boxplot(width=0.03, color="black", alpha=10) +
  scale_y_continuous(
    trans='log10',
    breaks = c(0, 1, 10, 100),
    guide = "axis_logticks",
    ) +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
          theme_classic()
dev.off()

# write master table
write.table(outtab_id3_TSS  , file=paste0(cellname, '.MergedFastq.node_itxTrans_chiadrop_nonPairwise.tsv'), sep='\t', quote=F, col.names=T, row.names=F)

