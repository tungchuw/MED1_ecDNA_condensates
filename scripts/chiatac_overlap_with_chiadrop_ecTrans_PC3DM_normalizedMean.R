#given 2 libraries, find differential iPET 
#trans itx: collect all anchors from libraries & replicates
#combine replicas

#merge the overlapping anchors into interaction nodes
#assign id for every node, thus the interaction between 2 nodes will have ID like node1_node2 
#Get a master table for all available samples iPET for each node1_node2 
#normalize with DEseq2
#find overlap with chiadrop anchors


#Rseurat
require(data.table)
require(ggplot2)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggrepel)
library(scales)


cellname='PC3DM+'

treatGroup='all'
# SE1Group='SE1'
# SE2Group='SE2'
# SE6Group='SE6'
# SE7Group='SE7'
# CombineGroup='combo_SE'


cdtsv=paste0('/net/nwgc/vol1/sharing/Wei_Lab/USERS/tungch/tjongh/ec/combine2Reps/',cellname,'.anchors_ecGEMs.annot.tsv')
chiadrop = read.delim(cdtsv)



# Directory paths
base_dir <- "./data/misc/"
cell_dir <- "./data/misc/PC3DM/"
script_dir <- "./scripts/"

# Input files
cdtsv <- paste0('./data/misc/PC3DM/chiadrop/', cellname, '.anchors_ecGEMs.annot.tsv')
sample_file <- file.path(cell_dir, paste0("chiatac/",cellname, "_all_itxlist.txt"))
samplePhenotype <- read.table(sample_file)
colnames(samplePhenotype) <- c("name","condition", "rep")
itxlibs <- samplePhenotype$name

celltreat <- 'combo_SE'
# cell_SE2 <- 'SE2'
# cell_SE6 <- 'SE6'
# cell_SE7 <- 'SE7'
# cell_combo_SE <- 'combo_SE'
# cell_combine_SE <- 'combine_SE'
cellcntrl <- 'vector'
compname <- paste0(celltreat, '_vs_', cellcntrl)
nsample = nrow(samplePhenotype)
comparedLibs  = split(samplePhenotype$name, samplePhenotype$condition)

bedtools='singularity run --bind /net/nwgc/vol1/nobackup/nocleanup/tungch/:/net/nwgc/vol1/nobackup/nocleanup/tungch/ /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools intersect -wo -a '
bedmerge='singularity run --bind /net/nwgc/vol1/nobackup/nocleanup/tungch/:/net/nwgc/vol1/nobackup/nocleanup/tungch/ /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools merge -i '
bedpetools='singularity run --bind /net/nwgc/vol1/nobackup/nocleanup/tungch/:/net/nwgc/vol1/nobackup/nocleanup/tungch/ /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools pairtopair -is -type both -a '

itxlibs = samplePhenotype$name;
itxsuf='.e500.clusters.trans.chiasig.gz'

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
itxDir = "/net/nwgc/vol1/nobackup/nocleanup/tungch/test/chiatac/ChiATAC_AT4-40/merged/" #modify accordingly
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

## get anchors
allitx = bind_rows(itxs) #combine all itx
anchors1 = allitx[,1:3]
anchors2 = allitx[,4:6]
colnames(anchors1) = colnames(anchors2) = c('CHR','Start','End')
anchors = rbind(anchors1, anchors2)
dim(anchors)
# [1] 25784480      3
rm(anchors1,anchors2)

#indexing the combined anchors
tmpbed = paste0(cellname, '.tmp.bed'); 
tmpout = paste0(cellname, '.tmp.out');
write.table(anchors, file=tmpbed, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0('sort -k1,1 -k2,2n ', tmpbed, ' | ', bedmerge, ' stdin > ', tmpout) 
system(com)
anchors = fread(tmpout, stringsAsFactors=F)
dim(anchors)
# [1] 951686      3
#colnames(anchors) =  c('CHR','Start','End'); anchors$ID = 1:nrow(anchors)
colnames(anchors) =  c('CHR','Start','End'); anchors$ID = paste0(anchors$CHR, ':', anchors$Start,'-',anchors$End)
anchorfile = paste0(cellname, '.anchors.bed')
write.table(anchors, file=anchorfile, sep='\t', quote=F, col.names=F, row.names=F)
com = paste0(bedtools, anchorfile, ' -b ', ecDNAfile, ' > ', tmpout)
system(com)
o = read.table(tmpout, stringsAsFactors=F)
dim(o)
# [1] 327   9
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

names(anchor.list)




# Merge to a matrix
alltab = merge(anchor.list[[1]],anchor.list[[2]], all=T)
alltab_2 = merge(alltab,anchor.list[[3]], all=T)

head(alltab_2)

library(tidyverse)
# Separate node_L and node_R
alltab_3 <- alltab_2 %>% 
  mutate(pair_L = sub("_.*", "", pairID),
         pair_R = sub(".*_", "", pairID)) %>%
  separate(pair_L, into = c("chr_L", "range_L"), sep = ":") %>% 
  separate(range_L, into = c("start_L", "end_L"), sep = "-") %>%
  separate(pair_R, into = c("chr_R", "range_R"), sep = ":") %>% 
  separate(range_R, into = c("start_R", "end_R"), sep = "-")

# filter chrM
table(alltab_3$chr_L != "chrM" & alltab_3$chr_R != "chrM")
  #  FALSE     TRUE 
  # 324514 11909108 
rawcounts = subset(alltab_3, chr_L != "chrM" & chr_R != "chrM", select = c(1:4))
dim(rawcounts)
# [1] 11909108       4

# Normalization
rawcounts[is.na(rawcounts)] = 0;
row.names(rawcounts) = rawcounts$pairID
rawcounts = rawcounts[,2:4]
dim(rawcounts)
# [1] 11909108        3


head(rawcounts)

normFactor = apply(rawcounts, 2, mean) #median wouldn't work since mock median = 0
normMat = t(apply(rawcounts, 1, function(x) x/normFactor)) 
colnames(normMat) = paste0('norm_', colnames(normMat))


DEtable = cbind(rawcounts, normMat); #check before merging that rownames are identical
DEtable$itxID = row.names(DEtable)
head(DEtable)




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
dim(mastertab)
# [1] 11909108       20
colnames(mastertab)
mastertab_filtered = mastertab


ec_ch = subset(mastertab_filtered, ecDNA == "L", select = c(7,1:6,10,13:20))
dim(ec_ch)
# [1] 25856    19
colnames(ec_ch)

# chr_ec
ch_ec = subset(mastertab_filtered, ecDNA == "R", select = c(7,1:6,11,13:20))
dim(ch_ec)
# [1] 313728     19

colnames(ch_ec)
head(ch_ec)

ch_ec$new_pairID = paste0(ch_ec$chr_R, ":" ,ch_ec$start_R,'-', ch_ec$end_R,"_",ch_ec$chr_L, ":" ,ch_ec$start_L,'-', ch_ec$end_L)
head(ch_ec)

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

dim(new.ch_ec)
head(new.ch_ec)
colnames(new.ch_ec)
#  [1] "itxID"                          "AT4.40.19m20.vector_iPET"      
#  [3] "AT4.40.1m2.SE1_iPET"            "AT4.40.9m10.combo_SE_iPET"     
#  [5] "norm_AT4.40.19m20.vector_iPET"  "norm_AT4.40.1m2.SE1_iPET"      
#  [7] "norm_AT4.40.9m10.combo_SE_iPET" "SErank_L"                      
#  [9] "TSS_L"                          "TSS_R"                         
# [11] "chr_L"                          "start_L"                       
# [13] "end_L"                          "chr_R"                         
# [15] "start_R"                        "end_R"  
colnames(ec_ch)
#  [1] "itxID"                          "AT4.40.19m20.vector_iPET"      
#  [3] "AT4.40.1m2.SE1_iPET"            "AT4.40.9m10.combo_SE_iPET"     
#  [5] "norm_AT4.40.19m20.vector_iPET"  "norm_AT4.40.1m2.SE1_iPET"      
#  [7] "norm_AT4.40.9m10.combo_SE_iPET" "SErank_L"                      
#  [9] "TSS_L"                          "TSS_R"                         
# [11] "chr_L"                          "start_L"                       
# [13] "end_L"                          "chr_R"                         
# [15] "start_R"                        "end_R" 

ecConnected.info = bind_rows(ec_ch, new.ch_ec)
dim(ecConnected.info)
# [1] 339584     16
head(ecConnected.info)


# save.image("temp.RData")



table(duplicated(ecConnected.info$itxID))
colnames(ecConnected.info)
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
dim(merged_ecConnected)
# [1] 339584     19

table(merged_ecConnected$TSS_R == ".")

# save.image(paste0('temp.RData'))
# # load(paste0(cellname, '.', compname,"_thcount",thcount, '.node_itxTrans.RData'))

colnames(merged_ecConnected)
#  [1] "itxID"                          "AT4.40.19m20.vector_iPET"      
#  [3] "AT4.40.1m2.SE1_iPET"            "AT4.40.9m10.combo_SE_iPET"     
#  [5] "norm_AT4.40.19m20.vector_iPET"  "norm_AT4.40.1m2.SE1_iPET"      
#  [7] "norm_AT4.40.9m10.combo_SE_iPET" "TSS_L"                         
#  [9] "TSS_R"                          "SErank_L"                      
# [11] "chr_L"                          "start_L"                       
# [13] "end_L"                          "chr_R"                         
# [15] "start_R"                        "end_R"     

dim(merged_ecConnected)
# [1] 339584     16

#Made table for chiadrop
ecConnected_regions <- merged_ecConnected[,c(14:16,11:13,1:10)]
head(ecConnected_regions)
dim(ecConnected_regions)
# [1] 339584     16


write.table(ecConnected_regions, file=tmpout, sep='\t', quote=F, col.names=F, row.names=F)
# Now tmpout contains chiapet ec-connected chromosomal anchors which will be picked if overlapped with chiadrop

#Overlap chiapet with chiadrop
chiadrop.P = subset(chiadrop
                    , ec == FALSE & PGI =='P'
                    , select=c(Chr, Start, End, nGEMs,MergeFragID, TSS)
                    )

colnames(chiadrop)

#counts
N.chiadrop <- nrow(chiadrop)
message("N.chiadrop: ", N.chiadrop)
# N.chiadrop: 14850
N.chiadrop.P <- nrow(chiadrop.P)
message("N.chiadrop.P: ", N.chiadrop.P)
#PC3.N.chiadrop.P: 10855

tmpout <- fread(tmpout, stringsAsFactors=F)
N.tmpout <- nrow(tmpout)
message("N.tmpout: ", N.tmpout)
#PC3.N.tmpout: 339584

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
#PC3_N.oo: 26105


colnames(oo) = c(colnames(ecConnected_regions), colnames(chiadrop.P), 'bpOverlap')

table(duplicated(oo$itxID))
# FALSE  TRUE 
# 18574  7531 

#reomve duplicated itxID
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
dim(as.data.frame(oo_unique))
# [1] 18574    23
outtab_id3_TSS <- as.data.frame(oo_unique)


for (j in 1:23) {
  count <- sum(grepl("\\b;\\b", outtab_id3_TSS[,j], perl = TRUE))
  if (count != 0) {
    print(paste0(colnames(outtab_id3_TSS)[j], ":", count))
  }
}



## unique gene number check
# separate
library(tidyr)
colnames(outtab_id3_TSS)
outtab_chiadrop_TSS <-  separate_rows(outtab_id3_TSS, chiadrop_TSS, sep = ";")
outtab_chiadrop_TSS = as.data.frame(outtab_chiadrop_TSS)
dim(outtab_chiadrop_TSS)
# [1] 39639    23
table(duplicated(outtab_chiadrop_TSS$chiadrop_TSS))
# FALSE  TRUE 
#  9212 30427  
length(unique(outtab_chiadrop_TSS$chiadrop_TSS))
# [1] 9212

outtab_TSS_R <-  separate_rows(outtab_id3_TSS, TSS_R, sep = ";")
outtab_TSS_R = as.data.frame(outtab_TSS_R)
dim(outtab_TSS_R)
# [1] 33184    23
table(duplicated(outtab_TSS_R$TSS_R))
# FALSE  TRUE 
#  9780 23404 
length(unique(outtab_TSS_R$TSS_R))
# [1] 9780



table(duplicated(outtab$chiadrop_TSS))
# FALSE  TRUE 
#  4190   669 
length(unique(outtab$chiadrop_TSS))
# [1] 4190



colnames(outtab_id3_TSS)
#  [1] "itxID"                          "Chr_R"                         
#  [3] "Start_R"                        "End_R"                         
#  [5] "Chr_L"                          "Start_L"                       
#  [7] "End_L"                          "AT4.40.19m20.vector_iPET"      
#  [9] "AT4.40.1m2.SE1_iPET"            "AT4.40.9m10.combo_SE_iPET"     
# [11] "norm_AT4.40.19m20.vector_iPET"  "norm_AT4.40.1m2.SE1_iPET"      
# [13] "norm_AT4.40.9m10.combo_SE_iPET" "TSS_L"                         
# [15] "TSS_R"                          "SErank_L"                      
# [17] "chiadrop_Chr"                   "chiadrop_Start"                
# [19] "chiadrop_End"                   "chiadrop_nGEMs"                
# [21] "chiadrop_MergeFragID"           "chiadrop_TSS"                  
# [23] "bpOverlap"
table(outtab_id3_TSS$SErank_L)
# using chiadrop genes
#     0     1     2    20    21     3    31    36     4     6    70 
# 14618  3771   971  2989   964  4859  1446   398  1249  2515   488 

# using betools 
#    0    1    2   20   21    3   31   36    4    6   70 
# 7978 2012  535 1630  514 2547  782  224  696 1386  270 

# crispr subset
outtab_id3_TSS_crispr = subset(outtab_id3_TSS, SErank_L %in% c(1,2,6,7))
dim(outtab_id3_TSS_crispr)
# [1] 3933   23
outtab_id3_TSS_SE1 = subset(outtab_id3_TSS, SErank_L == "1")
dim(outtab_id3_TSS_SE1)
# [1] 2012   23
# subset

Values_morm_vector_19n20 = as.numeric(outtab_id3_TSS_SE1$norm_AT4.40.19m20.vector_iPET[outtab_id3_TSS_SE1$norm_AT4.40.19m20.vector_iPET != 0])
Values_morm_SE1_1n2 = as.numeric(outtab_id3_TSS_SE1$norm_AT4.40.1m2.SE1_iPET[outtab_id3_TSS_SE1$norm_AT4.40.1m2.SE1_iPET != 0])
Values_morm_comboSE_9n10 = as.numeric(outtab_id3_TSS_SE1$norm_AT4.40.9m10.combo_SE_iPET[outtab_id3_TSS_SE1$norm_AT4.40.9m10.combo_SE_iPET != 0])

wilcox.test(Values_morm_SE1_1n2, 
            Values_morm_vector_19n20, 
            alternative='less')$p.value

# p-value = 1.369851e-158 SE1
wilcox.test(Values_morm_comboSE_9n10, 
            Values_morm_vector_19n20, 
            alternative='less')$p.value

# p-value = 1.640885e-175 SE1

nA_v = length(Values_morm_vector_19n20)

# [1] 578 SE1

nA_1 = length(Values_morm_SE1_1n2)

# [1] 783 SE1

nA_c = length(Values_morm_comboSE_9n10)
# [1] 8906
# [1] 1858 SE1,2,6,7
# [1] 904 SE1

cplabel = c(rep("vector", nA_v), rep("SE1", nA_1), rep("combo_SE", nA_c))
cp.df = data.frame(normFreq=c(Values_morm_vector_19n20, Values_morm_SE1_1n2, Values_morm_comboSE_9n10),  cond=factor(cplabel, levels=c("vector","SE1","combo_SE"))) #check if grpnames call is consistent

summary(Values_morm_vector_19n20)
summary(Values_morm_SE1_1n2)
summary(Values_morm_comboSE_9n10)
# > summary(Values_morm_vector_19n20)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   3.770   3.770   3.770   4.819   3.770 301.619 
# > summary(Values_morm_SE1_1n2)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.671   2.671   2.671   3.388   2.671 213.661 
# > summary(Values_morm_comboSE_9n10)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.316   2.316   2.316   2.918   2.316 317.247 

###########SE1,2,6,7
# > summary(Values_morm_vector_19n20)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   3.770   3.770   3.770   4.894   3.770 128.188 
# > summary(Values_morm_SE1_1n2)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.671   2.671   2.671   3.488   2.671  96.147 
# > summary(Values_morm_comboSE_9n10)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.316   2.316   2.316   3.017   2.316 101.890 

#########SE1
# > summary(Values_morm_vector_19n20)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   3.770   3.770   3.770   4.925   3.770 128.188 
# > summary(Values_morm_SE1_1n2)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.671   2.671   2.671   3.462   2.671  96.147 
# > summary(Values_morm_comboSE_9n10)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.316   2.316   2.316   2.984   2.316 101.890 

##Interaction violin
library(ggplot2)
library(ggbreak) 
library(scales)
library(gg.gap)
library(ggbeeswarm)

a=subset(cp.df)



pdf(paste0(cellname, '.MergedFastq.Trans.chiadrop_genes_nonPairwise_SE1filtered_bedtools.pdf'), width=6, height=6)
ggplot(a, aes(x=cond, y=normFreq, fill=cond)) + 
  geom_violin(width=0.7, bw =0.05) +
  scale_x_discrete(labels=c(paste0(unique(cplabel)[1]," (n = ",nA_v,")"), paste0(unique(cplabel)[2]," (n = ",nA_1,")"), paste0(unique(cplabel)[3]," (n = ",nA_c,")"))) +
  xlab("") +

  # geom_beeswarm() +
  # geom_quasirandom(dodge.width=.75, col=2) +
  # geom_point(position=position_jitterdodge(jitter.width = 0,dodge.width = 2)) +
  #scale_color_manual(values=c("#B099BF", "#8FBF88")) +
  scale_fill_manual(values=c("vector" = "#7fbf7b", "SE1" = "#af8dc3", "combo_SE" = "#af8dc3")) +
  geom_boxplot(width=0.03, color="black", alpha=10) +
  # geom_boxplot(width=0.15, color="black", alpha=0.2, outlier.shape = NA) +
  # ylim(0, 120) +
  scale_y_continuous(
    # labels = scientific,
    trans='log10',
    # trans = pseudo_log_trans(base = 10),
    breaks = c(0, 1, 10, 100),
    # minor_breaks = mb,
    guide = "axis_logticks",
    # labels = label_number(accuracy = 1),
    # minor_breaks = log10_minor_breaks()
    ) +
        #  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
        #       labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position="bottom",
        panel.spacing.x = unit(1.5, "lines")) +
          theme_classic()
dev.off()

# write master table
write.table(outtab_id3_TSS  , file=paste0(cellname, '.', compname, '.MergedFastq.node_itxTrans_chiadrop_nonPairwise.tsv'), sep='\t', quote=F, col.names=T, row.names=F)



