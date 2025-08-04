
library(circlize)
library(ComplexHeatmap)
library(gtools)
library(zoo)
library(gridBase)
library(dplyr)
library(stringr)
library(metafolio)
library(ggplot2)
require(gplots)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

bedmerge='singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools merge  -c 4 -o sum '
distiMerge='singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools merge -c 4 -o distinct '
bedsort='singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools sort -i'
bwprog='/net/nwgc/vol1/home/tungch/miniconda3/envs/bigwig/bin/bigWigSummary'
ncolor = 5
bordcol='#f5f5f5'

circosplot_EC_gene <- function(col.se='#9834eb', col.chip3='#f40000', col.chip2='#3acae0', col.chip1='#e0a33a', col.score.trans='#04cc65', col.score.cis='black', cex.se=0.7, geneplot=FALSE, maj.tick=1e9){
    circos.clear()
    circos.par("gap.degree"=1, start.degree = 90);
    om = circos.par("track.margin")
    oc = circos.par("cell.padding")
    chrnames = sapply(strsplit(as.vector(cytoband.df$V1),'_'),'[[',1);
    col3 = gg_color_hue(5)[-1]; names(col3) = c('chr7','chr8','chr12', 'chr1')
    circos.par(track.margin = c(om[1], 0), cell.padding = c(0, 0, 0, 0))
    circos.genomicInitialize(cytoband.df, tickLabelsStartFromZero=F, major.by=maj.tick, axis.labels.cex =1, labels.cex =0.6, sector.names=rep('',nrow(cytoband.df)))
    write.table(cytoband.df, file=paste0(RUN, '.cytoband.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")
    circos.track(ylim = c(0, .01), bg.col =  col3[chrnames], bg.border = NA, track.height = 0.02)
    circos.par(track.margin = c(om[1], 0), cell.padding = c(0, 0, 0, 0))
    genecolor = 1:nrow(genetable) %% ncolor + 1
    if(geneplot){
        circos.genomicLabels(genetable, labels.column = 4, side = "inside",col=genecolor, line_col=genecolor, cex=0.4, connection_height=convert_height(1, "mm"), padding=0.1)
        write.table(genetable, file=paste0(RUN, '.genetable.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")
    }else{
        circos.genomicLabels(genetable, labels.column = 4, side = "inside",col='white', line_col=genecolor, cex=0.1, connection_height=convert_height(2, "mm"), labels_height=0.05, padding=0.1)
    }

    #H3K27Ac Chipseq
    circos.genomicTrackPlotRegion(bwsig.df2[,c(1:4)], ylim = c(0, 17), stack = FALSE, panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, area =TRUE, border=NA,  type='l', lwd=2, col = col.chip3, ...)
    },  track.height = 0.12,  bg.border = bordcol)
        write.table(bwsig.df2, file=paste0(RUN, '.h3k27ac.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")

    if(ncol(bwsig.df2)>4 & geneplot==FALSE){
        circos.genomicTrackPlotRegion(bwsig.df2[,c(1:3,5)], ylim = c(0, 17), stack = FALSE, panel.fun = function(region, value, ...) {
            circos.genomicLines(region, value, area =TRUE, type='l', lwd=2, border=NA, col = col.chip1, ...)
        },  track.height = 0.12,  bg.border = NA)
    }

    #MED1 Chipseq
    circos.genomicTrackPlotRegion(bwsig.df[,c(1:4)], ylim = c(0, 7), stack = FALSE, panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, area =TRUE, border=NA,  type='l', lwd=2, col = col.chip2, ...)
    },  track.height = 0.12,  bg.border = bordcol)
        write.table(bwsig.df, file=paste0(RUN, '.med1.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")

    if(ncol(bwsig.df)>4 & geneplot==FALSE){
        circos.genomicTrackPlotRegion(bwsig.df[,c(1:3,5)], ylim = c(0, 7), stack = FALSE, panel.fun = function(region, value, ...) {
            circos.genomicLines(region, value, area =TRUE, type='l', lwd=2, border=NA, col = col.chip1, ...)
        },  track.height = 0.12,  bg.border = NA)
    }

    circos.genomicTrackPlotRegion(nodeScore.cis, ylim = c(0, 520), stack = FALSE, panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, area =FALSE, border=NA,  type='h', lwd=1, col = col.score.cis, ...)
    },  track.height = 0.15,  bg.border = bordcol)
    write.table(nodeScore.cis , file=paste0(RUN, '.cis_score.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")

    circos.genomicTrackPlotRegion(nodeScore.trans, ylim = c(0, 850), stack = FALSE, panel.fun = function(region, value, ...) {
       circos.genomicLines(region, value, area =FALSE, border=NA,  type='h', lwd=1, col = col.score.trans, ...)
    },  track.height = 0.12,  bg.border = bordcol)
    write.table(nodeScore.trans, file=paste0(RUN, '.trans_score.circos.tsv'), col.names=F, row.names=F, quote=F, sep="\t")

    message('Drawing links in circos')
    circos.genomicLink(bedL, bedR, col=ipet.col)
    circos.par(track.margin = om, cell.padding = oc)
    circos.clear()
}

#-------------------------------------------------------------------------------------------------------------------

intersectBed='singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools intersect -u -a '
intersectPair='singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools pairtobed -type both -a '
woIntersect='singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools intersect -wo -a'

RUN='ACD0004'
# RUN='ACD0005'
#RUN='ACD0013'
#RUN='ACD0014'
#RUN='ACD0015'
#RUN='ACD0016'
# RUN='ACD0034'
# RUN='ACD0035'

#replace .ecFragPair.txt with output of: 
chiadropDir='/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/chiadrop/'
chiapetDir='/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/ChIA-PIPE/'

run2names = c('PC3NCI60', 'PC3ATCC', 'COLO320DM','COLO320HSR', 'B168');#folder names /projects/wei-lab/proj-SuperEnhancer/ChIA-PIPE
run2names = c(rep(run2names[1:2], each=4), rep(run2names[3:5], each=2),  run2names[5], run2names[5], 'B171','B171')
names(run2names) =  paste0('ACD00',c('06','07',13:16,'11','12','04','05',17:20,'08', '23', 34:35)) #Chiadrop RUNs
itxFile <- "/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/ChIA-PIPE/HDtreat/ACP0016/ACP0016.DA.itx"

ecNames = c('COLO320DM', 'COLO320DM', 'PC3DM', 'PC3DM', 'B168','B171'); #prefix of files
names(ecNames) = c('COLO320DM', 'COLO320HSR', 'PC3ATCC', 'PC3NCI60', 'B168','B171'); #RUNs
seNames = c('COLO320DM', 'COLO320HSR', 'PC3DM-', 'PC3DM+', 'B168','B171');#filename .SErose.bed
names(seNames) = c('COLO320DM', 'COLO320HSR', 'PC3ATCC', 'PC3NCI60', 'B168','B171');#filename .SErose.bed

regiondir='/net/nwgc/vol1/nobackup/nocleanup/tungch/regions/'

medlibs = list()
medlibs[["PC3DM-"]] = c('A0035')
medlibs[["PC3DM+"]] = c('A0038')
medlibs[["COLO320DM"]] = c('A0026')
medlibs[["B168"]] = c('A0018')
medlibs[["B171"]] = c('A0041')
chipdir=paste0('/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/1x51_chipseq_processed/',med1.bwfiles[1],'/')
med1.bwfiles = medlibs[[seNames[run2names[RUN]]]]
bw1=paste0(chipdir,"fe_", med1.bwfiles[1], ".bw")
bw2=paste0("/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/1x51_chipseq_processed/A0010/normalized/fe_A0010.bw")
ecdnaRegfile = "B171.ecDNA.bed"
roseBed = paste0(seNames[run2names[RUN]] , '.SErose.bed')

#
genome='hg38'
inclCHR=paste0('chr',c(1:22,'X'))
allcyto = readRDS(system.file("extdata", "cytoband_list.rds", package = "circlize"))
cytoband.df = subset(allcyto[[genome]], V1 %in% inclCHR) 

#----------
options(scipen=999)
tmpbed = paste0(RUN,'.tmp.bed')
tmpout = paste0(RUN,'.tmpout.bed')
tmpcyto = paste0(RUN,'.tmpcyto.bed')
write.table(cytoband.df, file=tmpcyto, quote=F, sep="\t", col.names=F, row.names=F)

bpbin = 5000;
bpbin = 2500;
bpbin = 10000;
#Getting the interaction PET
bpbuff = 5000
ec.lib1 = read.table(ecdnaRegfile)[,1:3]
ec.lib1 = ec.lib1[order(ec.lib1$V1, ec.lib1$V2),]
message('Read ', ecdnaRegfile)
pti=95e6 #ticks

rangecmd = paste0('sort -k1,1 -k2,2n ', ecdnaRegfile, ' | singularity run /net/nwgc/vol1/nobackup/nocleanup/tungch/sifsh/chipseqtools.sif bedtools merge -d 1 -i stdin >', tmpout)
system(rangecmd)
ec.range <- read.table(tmpout, stringsAsFactors=F)
write.table(ec.range, file=tmpbed, quote=F, sep="\t", col.names=F, row.names=F)
bedcmd = paste(woIntersect, tmpbed, '-b', tmpcyto, '> ', tmpout)
system(bedcmd)
#REDFINE cytoband just ecDNA segments
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

#links for circos within ecDNA
bedcmd = paste0('cut -f1-7 ', itxFile, '|', intersectPair, 'stdin -b ', ecdnaRegfile, ' | cut -f1-7 |sort -k1,1 -k2,2n | uniq > ', tmpout)
system(bedcmd)
ec.itx = read.table(tmpout, stringsAsFactors=F)
message('Get ecDNA regions in ', itxFile )
nitx = nrow(ec.itx)
bedL = ec.itx[,1:3]
bedR = ec.itx[,4:6]
ipet = log10(ec.itx[,7])
message("Number of itx: ", nitx)
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

#Sum the overlapping fragments that come from different GEMs (fragments that made up interactions only)
load(paste0('/net/nwgc/vol1/sharing/Wei_Lab/USERS/tungch/tjongh/ec/test/anno_Rdata_COLO320DM_B171ecRegion/',RUN,'.annotChiaDrop.ecPol2.RData'))
message('Loading fragments from ', paste0('/net/nwgc/vol1/sharing/Wei_Lab/USERS/tungch/tjongh/ec/test/anno_Rdata_COLO320DM_B171ecRegion/',RUN,'.annotChiaDrop.ecPol2.RData'))
ecfrag = subset(fragdf, ecFrag==TRUE, select=c(GEMid, frag, GEMtype)) #fragdf is already non-singleton with frag rs3
eclist = split(as.vector(ecfrag$frag), as.vector(ecfrag$GEMid))
ecfrag.str = as.vector(unlist(eclist))
ecfrag = gsub('-','\t',gsub(':', '\t',ecfrag.str))
itxBed=paste(ecfrag, "1", sep='\t')
write.table(itxBed, file=tmpbed, quote=F, col.names=F, row.names=F)
bedcmd = paste(bedsort, tmpbed, '|', bedmerge, '|', woIntersect, 'stdin -b ', tmpcyto,  ' >', tmpout)
system(bedcmd)
itxnodes= read.table(tmpout, stringsAsFactors=F)
d = itxnodes$V3-itxnodes$V2
nodeScore.cis = data.frame(Chr=itxnodes$V8, Start=itxnodes$V2, End=itxnodes$V3, Score=itxnodes$V4/d*1000) # can be written as a table

#chromosomal partner score
ectrans = subset(fragdf, GEMtype=='ecTrans', select=c(GEMid, frag, ecFrag))
ectrans.ec = subset(ectrans, ecFrag==TRUE)
s = strsplit(as.vector(ectrans.ec$GEMid),'-')
transBC = sapply(s, '[[',4); #barcode of trans ec
ectrans.ec.bed = gsub('-','\t',gsub(':', '\t',ectrans.ec$frag)) #bed format of ec fragments that do trans
ecChr = subset(ectrans, ecFrag==FALSE ) #chromosomal fragments
s = strsplit(as.vector(ecChr$GEMid), '-')
chrBC =  sapply(s, '[[',4); #barcode of chromosomal in trans ec 
ecChrlist = split(as.vector(ecChr$frag), chrBC)
ecChr.nfrag = sapply(ecChrlist, length) #nfrag of chromosomal per GEM
nChrfrag_transGEM = ecChr.nfrag[transBC]
ectrans.ec.nfrag.bed = paste(ectrans.ec.bed, nChrfrag_transGEM, sep='\t')
write.table(ectrans.ec.nfrag.bed, file=tmpbed, quote=F, col.names=F, row.names=F)
bedcmd = paste(bedsort, tmpbed, '|', bedmerge, '|', woIntersect, 'stdin -b ', tmpcyto,  ' >', tmpout)
system(bedcmd)
itxnodes= read.table(tmpout, stringsAsFactors=F)
d = itxnodes$V3-itxnodes$V2 #for normalization
nodeScore.trans = data.frame(Chr=itxnodes$V8, Start=itxnodes$V2, End=itxnodes$V3, Score=itxnodes$V4/d*1000) # can be written as a table

#colors
transpare = 0.7
ipet.bar = c(0.47, 1, 1.5)
col_fun = colorRamp2(ipet.bar, c("#daf542", "#f5aa42", "#e30e0e"), transparency=transpare)
ipet.col = col_fun(ipet);
lgd_link = Legend(at=ipet.bar, col_fun = col_fun, title_position = "topcenter", title = "log10(iPET count)", direction = "horizontal", border = "black");
lgd_list_horizontal = packLegend(lgd_link);
col_fun_plot = colorRamp2(ipet.bar, c("#daf542", "#f5aa42", "#e30e0e"), transparency=0)
lgd_link_plot = Legend(at=ipet.bar, col_fun = col_fun_plot, title_position = "topcenter", title = "log10(iPET count)", direction = "horizontal", border = "black");
lgd_list_horizontal_plot = packLegend(lgd_link_plot);

#call chipseq on bins, to reduce the difficulty plotting
#bin it
bwsig.df = NULL;
for (i in 1:nrow(ec.range)){
    ch = as.vector(ec.range[i,1]);
    p1=ec.range[i,2];
    p2=ec.range[i,3];
    nbin = as.numeric(floor((p2-p1)/bpbin)); #number of bins in an ec range
    start.pos = seq(p1, by=bpbin, length.out=nbin)
    stop.pos = start.pos + bpbin-1
    last.start = stop.pos[nbin]+1 #for excess bin less than bpbin
    last.stop = p2  #for excess bin less than bpbin
    exd = last.stop - last.start
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
    nbin = as.numeric(floor((p2-p1)/bpbin)); #number of bins in an ec range
    start.pos = seq(p1, by=bpbin, length.out=nbin)
    stop.pos = start.pos + bpbin-1
    last.start = stop.pos[nbin]+1 #for excess bin less than bpbin
    last.stop = p2  #for excess bin less than bpbin
    exd = last.stop - last.start
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
i=grep('^MIR[0-9]',genetable$Gene)
if (length(i) > 0 ){genetable = genetable[-i,]}
system(paste('rm', tmpbed, tmpcyto, tmpout))

#Circos
circle_size = unit(1, "snpc") # snpc unit gives you a square region
#Plotting
vp <- viewport(width=0.95, height=0.95)
plot.new()
par(cex=1.3)
bplab=paste0(bpbin/1000,'kb')

pdf(paste0(RUN, '.', seNames[run2names[RUN]], ".ecCircos.",bplab,"_",med1.bwfiles[1],"_with_ChIAPET_trans", transpare,".pdf"), height=8, width=8, pointsize=15);
circosplot_EC_gene(geneplot=F, maj.tick=1E10, col.chip1='#b00448')
pushViewport(viewport(x=0.8, y=0.05, width=0.4, height=0.05, name="Leg", just=c("center","bottom")))
grid.draw(lgd_list_horizontal_plot);
circosplot_EC_gene(geneplot=T, maj.tick=1E10, col.chip1='#b00448')
pushViewport(viewport(x=0.7, y=0.19, width=0.4, height=0.05, name="Leg", just=c("center","bottom")))
grid.draw(lgd_list_horizontal_plot);
upViewport(0);
dev.off()
