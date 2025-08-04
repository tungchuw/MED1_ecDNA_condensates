#Draw ecDNA region ONLY circos with tracks:
#1. Fragment counts chiadrop ec-chr barplots (rs3, fragdf type==ecTrans in proj-SuperEnhancer/chiadrop/*annotChiaDrop.ecPol2.RData)
#   a=subset(fragdf, GEMtype=='ecTrans' & ecFrag==FALSE, select=c(frag, Nread))
#2. sum fragments ec
#3. MED1 profile
#4. Super enhancer annotation
#5. Gene names (just ticks and span)

library(circlize)
library(ComplexHeatmap)
library(gtools)
library(zoo)
library(gridBase)
library(dplyr)
library(stringr)
require(gplots)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

bedtools='singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools '
bedmerge='singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools merge  -c 4 -o sum '
distiMerge='singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools merge -c 4 -o distinct '
bedsort='singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools sort -i'

bwprog='/net/nwgc/vol1/home/tungch/miniconda3/envs/bigwig/bin/bigWigSummary'

ncolor = 5
bordcol='#f5f5f5'

#col.se='#9834eb'; col.chip2='#3acae0'; col.chip1='#e0a33a'; col.score.cis='#04cc65'; col.score.trans='black'; cex.se=0.6; geneplot=FALSE; maj.tick=5e6
circosplot_EC_gene <- function(ynum, col.se='#9834eb', col.chip3='#f40000', col.chip2='#3acae0', col.chip1='#e0a33a', col.score.trans='#04cc65', col.score.cis='black', cex.se=0.7, geneplot=FALSE, maj.tick=1e9, trh=0.1){

    circos.clear()

    circos.par("gap.degree"=1, start.degree = 90);
    om = circos.par("track.margin")
    oc = circos.par("cell.padding")

    chrnames = sapply(strsplit(as.vector(cytoband.df$V1),'_'),'[[',1);
    col3 = gg_color_hue(4)[-1]; names(col3) = c('chr7','chr8','chr12')
    circos.par(track.margin = c(om[1], 0), cell.padding = c(0, 0, 0, 0))

    #circos.genomicInitialize(cytoband.df, tickLabelsStartFromZero=F, major.by=maj.tick, axis.labels.cex =1, labels.cex =1.2)
    circos.genomicInitialize(cytoband.df, tickLabelsStartFromZero=F, major.by=maj.tick, axis.labels.cex =1, labels.cex =0.6, sector.names=rep('',nrow(cytoband.df)))
    write.table(cytoband.df, file=paste0(cellname, '.cytoband.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")
    circos.track(ylim = c(0, .01), bg.col =  col3[chrnames], bg.border = NA, track.height = 0.02)
    circos.par(track.margin = c(om[1], 0), cell.padding = c(0, 0, 0, 0))
    genecolor = 1:nrow(genetable) %% ncolor + 1
    if(geneplot){
        #i=which(sapply(as.vector(genetable$Gene), function(x){ nchar(x)}) < 5)
        #genetable$Label[i] = as.vector(genetable$Gene[i])
        #genetable = genetable[order(genetable$Label, decreasing=T),]
        circos.genomicLabels(genetable, labels.column = 4, side = "inside",col=genecolor, line_col=genecolor, cex=0.4, connection_height=convert_height(1, "mm"), padding=0.1)
        write.table(genetable, file=paste0(cellname, '.genetable.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")
    }else{
        circos.genomicLabels(genetable, labels.column = 4, side = "inside",col='white', line_col=genecolor, cex=0.1, connection_height=convert_height(2, "mm"), labels_height=0.05, padding=0.1)
    }

    #H3K27Ac Chipseq
    circos.genomicTrackPlotRegion(bwsig.df2[,c(1:4)], ylim = c(0, 12), stack = FALSE, panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, area =TRUE, border=NA,  type='l', lwd=2, col = col.chip3, ...)
    },  track.height = 0.12,  bg.border = bordcol)
        write.table(bwsig.df2, file=paste0(cellname, '.h3k27ac.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")

    if(ncol(bwsig.df2)>4 & geneplot==FALSE){
        circos.genomicTrackPlotRegion(bwsig.df2[,c(1:3,5)], ylim = c(0, 12), stack = FALSE, panel.fun = function(region, value, ...) {
            circos.genomicLines(region, value, area =TRUE, type='l', lwd=2, border=NA, col = col.chip1, ...)
        },  track.height = 0.12,  bg.border = NA)
    }


    #MED1 Chipseq
    circos.genomicTrackPlotRegion(bwsig.df[,c(1:4)], ylim = c(0, 30), stack = FALSE, panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, area =TRUE, border=NA,  type='l', lwd=2, col = col.chip2, ...)
    },  track.height = 0.12,  bg.border = bordcol)
        write.table(bwsig.df, file=paste0(cellname, '.med1.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")

    if(ncol(bwsig.df)>4 & geneplot==FALSE){
        circos.genomicTrackPlotRegion(bwsig.df[,c(1:3,5)], ylim = c(0, 30), stack = FALSE, panel.fun = function(region, value, ...) {
            circos.genomicLines(region, value, area =TRUE, type='l', lwd=2, border=NA, col = col.chip1, ...)
        },  track.height = 0.12,  bg.border = NA)
    }

    #SuperEnhancer ranks
    circos.genomicLabels(seReg[,c(1:3,5)] , labels.column = 4, side = "inside",col=col.se, line_col = 'black', cex=cex.se, connection_height=convert_height(1, "mm"), labels_height=0.06 )

    #circos.genomicTrack(ec.lib1, ylim=c(0,0.5),
    #    panel.fun = function(region, value, ...) {
    #    circos.genomicRect(region, value, col="red", border='red')
    #}, track.height = 0.02, bg.border = 'grey')

    # circos.genomicTrackPlotRegion(ec.mergeFrag, stack = FALSE, panel.fun = function(region, value, ...) {
    #     circos.genomicLines(region, value, area =FALSE, border=NA,  type='h', lwd=1, col = col.chip1, ...)
    # },  track.height = trh ,  bg.border = bordcol)
    # write.table(ec.mergeFrag , file=paste0(cellname, '.ec.mergeFrag.nGEMs.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")
    nodeScore.trans_filtered <- nodeScore.trans[nodeScore.trans$Score > ynum, ]

    circos.genomicTrackPlotRegion(nodeScore.trans_filtered, ylim = c(0, 250), stack = FALSE, panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, area =FALSE, border=NA,  type='h', lwd=1, col = col.score.trans, ...)
    },  track.height = trh ,  bg.border = bordcol)
    write.table(nodeScore.trans_filtered, file=paste0(cellname, '.trans_score_', ynum,'.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")

    # circos.genomicTrackPlotRegion(nodeScore.cis, stack = FALSE, panel.fun = function(region, value, ...) {
    #     circos.genomicLines(region, value, area =FALSE, border=NA,  type='h', lwd=1, col = col.score.cis, ...)
    # },  track.height = 0.15,  bg.border = bordcol)
    # write.table(nodeScore.cis , file=paste0(cellname, '.cis_score.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")

    #Links
    if(!geneplot){
        message('Drawing links in circos')
        circos.genomicLink(bedL, bedR, col=ipet.col);
    }

    circos.par(track.margin = om, cell.padding = oc)
    circos.clear()
}

#-------------------------------------------------------------------------------------------------------------------

intersectBed='singularity run --bind /net/nwgc/vol1/nobackup/nocleanup/tungch/:/net/nwgc/vol1/nobackup/nocleanup/tungch/ /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools intersect -u -a '
intersectPair='singularity run --bind /net/nwgc/vol1/nobackup/nocleanup/tungch/:/net/nwgc/vol1/nobackup/nocleanup/tungch/ /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools pairtobed -type both -a '
woIntersect='singularity run --bind /net/nwgc/vol1/nobackup/nocleanup/tungch/:/net/nwgc/vol1/nobackup/nocleanup/tungch/ /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools intersect -wo -a'

expname = 'AT4-36-m4'
cellname = 'PC3DM+'
regiondir='/net/nwgc/vol1/nobackup/nocleanup/tungch/regions/'
chipdir='/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/1x51_chipseq_processed/'
dataDir='/net/nwgc/vol1/nobackup/nocleanup/tungch/test/chiatac/AT4-36-merge/AT4-36-m4/ecTrans_DA/'
# chiadropDir='/net/nwgc/vol1/sharing/Wei_Lab/USERS/tungch/tjongh/ec/combine2Reps/'
itxFile <- paste0(dataDir, 'PC3DM+.chiatac.AT4-36-m4.annot_noM.tsv')#

#Files:
#peakfile=paste0(regiondir, cellname , ".Pol2_chipseq.narrowPeak")

medlibs = list()
medlibs[["PC3DM-"]] = c('A0035')
medlibs[["PC3DM+"]] = c('A0038')
medlibs[["COLO320DM"]] = c('A0026')
medlibs[["B168"]] = c('A0018')
med1.bwfiles = medlibs[[cellname]]
bw1=paste0(chipdir, med1.bwfiles[1],"/fe_", med1.bwfiles[1], ".bw")
bw2=paste0("/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/1x51_chipseq_processed/A0015/normalized/fe_A0015.bw")
bw1
bw2
ecdnaRegfile = paste0(cellname , '.ecDNA.bed')
roseBed = paste0(cellname , '.SErose.bed')

#
genome='hg38'
inclCHR=paste0('chr',c(1:22,'X'))
#allcyto = readRDS('/pod/2/wei-lab/USERS/tjongh/sumner/conda/envs/r-env/lib/R/library/circlize/extdata/cytoband_list.rds') #this contains hg38 as well
allcyto = readRDS(system.file("extdata", "cytoband_list.rds", package = "circlize"))
cytoband.df = subset(allcyto[[genome]], V1 %in% inclCHR) 


#----------
options(scipen=999)
tmpbed = paste0(cellname,'.tmp.bed')
tmpout = paste0(cellname,'.tmpout.bed')
tmpcyto = paste0(cellname,'.tmpcyto.bed')

write.table(cytoband.df, file=tmpcyto, quote=F, sep="\t", col.names=F, row.names=F)

bpbin = 10000;
#Getting the interaction PET
bpbuff = 5000
ec.lib1 = read.table(ecdnaRegfile)[,1:3]
ec.lib1 = ec.lib1[order(ec.lib1$V1, ec.lib1$V2),]
message('Read ', ecdnaRegfile)

pti=95e6 #ticks


rangecmd = paste0('sort -k1,1 -k2,2n ', ecdnaRegfile, '|' ,bedtools, 'merge -d 1 -i stdin >', tmpout)
system(rangecmd)
ec.range <- read.table(tmpout, stringsAsFactors=F)
write.table(ec.range, file=tmpbed, quote=F, sep="\t", col.names=F, row.names=F)

#Create cytoband inside ec
bedcmd = paste(woIntersect, tmpbed, '-b', tmpcyto, '> ', tmpout)
system(bedcmd)
cytoband.df = read.table(tmpout, stringsAsFactors=F)
cytoband.df = cytoband.df[,c(1:3,7,8)]
pad0 = str_pad(1:nrow(cytoband.df), width=2, pad="0")
cytochr = paste0(cytoband.df$V1,'_',pad0)
write.table(cbind(cytoband.df[,1:3],cytochr), file=tmpcyto, quote=F, sep="\t", col.names=F, row.names=F)
cytoband.df = data.frame(V1=cytochr, V2=cytoband.df[,2]-bpbuff, V3=cytoband.df[,3]+bpbuff)
#cytoband has extra bpbuff
nband = nrow(cytoband.df)
gaps = cytoband.df$V2[2:nband] - cytoband.df$V3[1:(nband-1)]
igaps =  which(gaps < 1)
if (length(igaps > 0)){message("Cytoband gap distance < 1 : ",igaps )}


#Gene table
geneSource='hg38.Genes.itxAnnotation.bed'
gtext="| grep -E 'TES|TSS'"
awktxt="| awk -F':' '{print $1,$2}'"
bedcmd = paste(intersectBed, geneSource,' -b', ecdnaRegfile, gtext, awktxt, '>', tmpout)
system(bedcmd)

gtmp = read.table(tmpout, stringsAsFactors=F)
gnames = unique(gtmp$V5)
ng = length(gnames)
ach = character(ng)
amin = numeric(ng)
amax = numeric(ng)
names(ach) = names(amin) = names(amax) = gnames
for ( g in gnames ){
    a = subset(gtmp, V5==g)
    ach[g] = unique(a$V1)
    amin[g] = min(a$V2)
    amax[g] = max(a$V3)
}

genetable = data.frame(Chr=ach, Start=amin,End=amax,Gene=gnames ) 
genetable0 = genetable
write.table(genetable, file=tmpbed, quote=F, sep="\t", col.names=F, row.names=F)
bedcmd = paste(woIntersect, tmpbed, '-b', tmpcyto, '> ', tmpout)
system(bedcmd)
g = read.table(tmpout, stringsAsFactors=F)
genetable = data.frame(Chr=g$V8, Start=g$V2,End=g$V3,Gene=g$V4)

#------

#Super Enhancers
bedcmd = paste(woIntersect, roseBed,' -b', tmpcyto,'>', tmpout)
system(bedcmd)
g = read.table(tmpout, stringsAsFactors=F)
message('Get ecDNA regions in ', roseBed )
seReg = data.frame(Chr=g$V8, Start=g$V2,End=g$V3, Rank=g$V4)
seReg$value = (1/seReg$Rank)/sum(1/seReg$Rank)
seReg=seReg[,c(1:3,5,4)]

#links for circos within ecDNA
#bedcmd = paste0('cut -f1-7 ', itxFile, '|', intersectPair, 'stdin -b ', ecdnaRegfile, ' | cut -f1-7 |sort -k1,1 -k2,2n | uniq > ', tmpout)
#system(bedcmd)
oitx = read.table(itxFile, header=T, sep="\t", comment.char='', check.names=FALSE, stringsAsFactors=FALSE);
message('Get ecDNA regions in ', itxFile )
ec.itx = subset(oitx, ecDNA == 'LR', select=1:7) #ec-ec
nitx = nrow(ec.itx) #for the circos links
message("Number of ec-ec itx: ", nitx)
ec.itx = ec.itx[order(ec.itx$iPET),]
ecL  = subset(oitx, ecDNA == 'L', select=1:7) #ec-chr
ecR  = subset(oitx, ecDNA == 'R', select=1:7) #chr-ec
midL = (ec.itx$startL + ec.itx$endL)/2
midR = (ec.itx$startR + ec.itx$endR)/2
dec = midR-midL
dmin = 100000 #drawing links above dmin
isel = which(dec > dmin);
#[,c(4:6,1:3,7)] #ec-chr

bedL = ec.itx[isel,1:3]
bedR = ec.itx[isel,4:6]
ipet = log2(ec.itx[isel,7])
#rename chr anchors:
write.table(bedL, file=tmpbed, quote=F, sep="\t", col.names=F, row.names=F)
bedcmd = paste(woIntersect, tmpbed,' -b', tmpcyto,'>', tmpout)
system(bedcmd)
g = read.table(tmpout, stringsAsFactors=F)
bedL = g[,c(7,2:3)]
write.table(bedR, file=tmpbed, quote=F, sep="\t", col.names=F, row.names=F)
bedcmd = paste(woIntersect, tmpbed,' -b', tmpcyto,'>', tmpout)
system(bedcmd)
g = read.table(tmpout, stringsAsFactors=F)
bedR = g[,c(7,2:3)]
#colnames(bedR) = colnames(bedL) #already the same

#Get the chiatac itx freq
o = ec.itx; 
colnames(o) = colnames(ecL) = colnames(ecR) = c('Chr','Start','End','Chr','Start','End','ipet')
chiatac.cis = rbind(o[,c(1:3,7)],o[,4:7]) #ec parts
chiatac.trans = rbind(ecL[,c(1:3,7)], ecR[,4:7]) #ec parts
write.table(chiatac.cis, file=tmpbed, quote=F, col.names=F, row.names=F, sep='\t')
bedcmd = paste(bedsort, tmpbed, '|', bedmerge, '|', woIntersect, 'stdin -b ', tmpcyto,  ' >', tmpout)
system(bedcmd)
itxnodes= read.table(tmpout, stringsAsFactors=F)
d = itxnodes$V3-itxnodes$V2
#The score is defined as number of piled up tags per kb
nodeScore.cis = data.frame(Chr=itxnodes$V8, Start=itxnodes$V2, End=itxnodes$V3, Score=itxnodes$V4/d*1000) # can be written as a table
#chromosomal partner score
write.table(chiatac.trans, file=tmpbed, quote=F, col.names=F, row.names=F, sep='\t')
bedcmd = paste(bedsort, tmpbed, '|', bedmerge, '|', woIntersect, 'stdin -b ', tmpcyto,  ' >', tmpout)
system(bedcmd)
itxnodes= read.table(tmpout, stringsAsFactors=F)
d = itxnodes$V3-itxnodes$V2
nodeScore.trans = data.frame(Chr=itxnodes$V8, Start=itxnodes$V2, End=itxnodes$V3, Score=itxnodes$V4/d*1000) # can be written as a table


# #Get Chiadrop interaction freq
# mergeFrag = read.delim(paste0(chiadropDir, cellname, '.anchors_ecGEMs.annot.tsv'), stringsAsFactors=F)
# ec.mergeFrag = subset(mergeFrag, ec==TRUE, select=c(Chr,Start,End,nGEMs))
# write.table(ec.mergeFrag, file=tmpbed, quote=F, sep="\t", col.names=F, row.names=F)
# bedcmd = paste(woIntersect, tmpbed,' -b', tmpcyto,'>', tmpout)
# system(bedcmd)
# g = read.table(tmpout, stringsAsFactors=F)
# ec.mergeFrag = g[,c(8,2:4)]


#colors
#ipet.bar = c(0.47, 1, 1.5) #for log10
#ipet.bar = c(2, 4, 6) #for log2
# ipet.bar = c(2, 4, 6, 8)
# col_fun = colorRamp2(ipet.bar, c("#E7E7E7", "#daf542", "#f5aa42", "#e30e0e"), transparency=0.7)
ipet.bar = c(2, 3, 5, 6)
col_fun = colorRamp2(ipet.bar, c("#E7E7E7", "#daf542", "#f5aa42", "#e30e0e"), transparency=0.7)
ipet.col = col_fun(ipet);
lgd_link = Legend(at=ipet.bar, col_fun = col_fun, title_position = "topcenter", title = "log2(PET count)", direction = "horizontal", border = "black");
lgd_list_horizontal = packLegend(lgd_link);


# Just for ploting scale bar
col_fun_plot = colorRamp2(ipet.bar, c("#E7E7E7", "#daf542", "#f5aa42", "#e30e0e"), transparency=0)
lgd_link_plot = Legend(at=ipet.bar, col_fun = col_fun_plot, title_position = "topcenter", title = "log10(iPET count)", direction = "horizontal", border = "black");
lgd_list_horizontal_plot = packLegend(lgd_link_plot);



#call chipseq on bins, to reduce the difficulty plotting
#bin it
bwsig.df = NULL;
for (i in 1:nrow(ec.range)){
    ch = as.vector(ec.range[i,1]);
    p1=ec.range[i,2];
    p2=ec.range[i,3];
    #tagname = paste0(ch,'_',p1)
    nbin = as.numeric(floor((p2-p1)/bpbin)); #number of bins in an ec range
    start.pos = seq(p1, by=bpbin, length.out=nbin)
    stop.pos = start.pos + bpbin-1
    last.start = stop.pos[nbin]+1 #for excess bin less than bpbin
    last.stop = p2  #for excess bin less than bpbin
    exd = last.stop - last.start
    #call bigwig
    message('Reading signal: ', bw1, ' and excess bin size: ', exd)
    chipcom = paste(bwprog, bw1, ch, p1, stop.pos[nbin], nbin)
    sig1 <-  as.numeric(unlist(strsplit(system(chipcom, intern=T), '\t')));
    sig1[is.na(sig1)] = 0
    chipcom = paste(bwprog, bw1, ch, last.start, last.stop, 1)
    if(exd > 10){
        sig2 <- as.numeric(system(chipcom, intern=T))
        sig2[is.na(sig2)] = 0
        bwsig.df = bind_rows(bwsig.df, data.frame(Chr=ch,Start=c(start.pos, last.start), End=c(stop.pos, last.stop), Med1=c(sig1,sig2)))
    }else{
        bwsig.df = bind_rows(bwsig.df, data.frame(Chr=ch,Start=start.pos, End=stop.pos, Med1=sig1))
    }
}
write.table(bwsig.df, file=tmpbed, quote=F, sep="\t", col.names=F, row.names=F)
bedcmd = paste(woIntersect, tmpbed,' -b', tmpcyto,'>', tmpout)
system(bedcmd)
g = read.table(tmpout, stringsAsFactors=F)
bwsig.df  = g[,c(8,2:4)]


#call chipseq on bins, to reduce the difficulty plotting
#bin it
bwsig.df2 = NULL;
for (i in 1:nrow(ec.range)){
    ch = as.vector(ec.range[i,1]);
    p1=ec.range[i,2];
    p2=ec.range[i,3];
    #tagname = paste0(ch,'_',p1)
    nbin = as.numeric(floor((p2-p1)/bpbin)); #number of bins in an ec range
    start.pos = seq(p1, by=bpbin, length.out=nbin)
    stop.pos = start.pos + bpbin-1
    last.start = stop.pos[nbin]+1 #for excess bin less than bpbin
    last.stop = p2  #for excess bin less than bpbin
    exd = last.stop - last.start
    #call bigwig
    message('Reading signal: ', bw2, ' and excess bin size: ', exd)
    chipcom = paste(bwprog, bw2, ch, p1, stop.pos[nbin], nbin)
    sig1 <-  as.numeric(unlist(strsplit(system(chipcom, intern=T), '\t')));
    sig1[is.na(sig1)] = 0
    chipcom = paste(bwprog, bw2, ch, last.start, last.stop, 1)
    if(exd > 10){
        sig2 <- as.numeric(system(chipcom, intern=T))
        sig2[is.na(sig2)] = 0
        bwsig.df2 = bind_rows(bwsig.df2, data.frame(Chr=ch,Start=c(start.pos, last.start), End=c(stop.pos, last.stop), Med1=c(sig1,sig2)))
    }else{
        bwsig.df2 = bind_rows(bwsig.df2, data.frame(Chr=ch,Start=start.pos, End=stop.pos, Med1=sig1))
    }
}
write.table(bwsig.df2, file=tmpbed, quote=F, sep="\t", col.names=F, row.names=F)
bedcmd = paste(woIntersect, tmpbed,' -b', tmpcyto,'>', tmpout)
system(bedcmd)
g = read.table(tmpout, stringsAsFactors=F)
bwsig.df2  = g[,c(8,2:4)]






#reduce genetable
#a=subset(bwsig.df, Med1>quantile(bwsig.df$Med1,prob=0.8))
#colnames(a)[1:3] = colnames(bedL)
#a = rbind(a[,1:3], itxBed)
#write.table(a, file=tmpbed, quote=F, sep="\t", col.names=F, row.names=F)
i=grep('^MIR[0-9]',genetable$Gene)

if (length(i) > 0 ){genetable = genetable[-i,]}

#system(paste('rm', tmpbed, tmpcyto, tmpout))

#Circos
circle_size = unit(1, "snpc") # snpc unit gives you a square region

#Plotting
vp <- viewport(width=0.95, height=0.95)
plot.new()
par(cex=1.3)
#grid.newpage()
#png(paste0(lib1, "_ecCircos.labelGenes.png"), height=600, width=600, pointsize=15);
bplab=paste0(bpbin/1000,'kb')


library(metafolio)
library(ggplot2)

ynum = 20

pdf(paste0(cellname,'_',expname, ".chiatac_ecCircos.",bplab, "_yMin_",ynum,"_scaled_noM.pdf"), height=8, width=8, pointsize=15);
circosplot_EC_gene(ynum = ynum, geneplot=F, maj.tick=1E10, col.chip1='#b00448')
pushViewport(viewport(x=0.8, y=0.05, width=0.4, height=0.05, name="Leg", just=c("center","bottom")))
grid.draw(lgd_list_horizontal_plot);
circosplot_EC_gene(ynum = ynum, geneplot=T, maj.tick=1E10, col.chip1='#b00448')
#
# circosplot_EC_gene(geneplot=F, maj.tick=1E10, col.chip1='#b00448')
# circosplot_EC_gene(geneplot=T, maj.tick=1E10, col.chip1='#b00448')
# #
##Put legend
pushViewport(viewport(x=0.8, y=0.05, width=0.4, height=0.05, name="Leg", just=c("center","bottom")))
grid.rect(gp=gpar(col="blue")) #checking only
grid.draw(lgd_list_horizontal_plot);
upViewport(0);

dev.off()

#source('ecDNAonly_circlize.R')




