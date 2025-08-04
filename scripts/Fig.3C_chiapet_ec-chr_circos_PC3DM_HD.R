#Draw ecDNA region ONLY circos with tracks:

library(circlize)
library(ComplexHeatmap)
library(gtools)
library(zoo)
library(gridBase)
library(dplyr)
library(stringr)
require(gplots)

genome='hg38'
inclCHR=paste0('chr',c(1:22,'X'))
allcyto = readRDS(system.file("extdata", "cytoband_list.rds", package = "circlize"))
cytoband.df = subset(allcyto[[genome]], V1 %in% inclCHR) 

# function
process_samples <- function(itxs,cell,group,min,max,max_color) {
    print(paste0(group,": ",nrow(itxs)))
    bedL = itxs[,1:3]
    bedR = itxs[,4:6]
    ipet = itxs[,7]
    breaks <- c(0,1:10,max(itxs[,7])+1)
    interval <- cut(itxs[,7], breaks = breaks, right = FALSE)
    frequency_table <- table(interval)
    print(frequency_table)
    # line color
    ipet.bar = c(min,2.75,max)
    col_fun = colorRamp2(ipet.bar, c("#E7E7E7","orange", max_color), transparency= 0.5)
    ipet.col = col_fun(ipet);
    # Scale bar
    lgd_link = Legend(at=ipet.bar, col_fun = col_fun, title_position = "topcenter", title = "iPET count", direction = "horizontal");
    lgd_list_horizontal = packLegend(lgd_link);
    # plot
    pdf(paste0(cell,"_Chia-PET_ecSE1_",group,"_min_",min,"_max_",max,"_","_scale_color_",max_color,"_updated.pdf"), height=8, width=8, pointsize=15);
        circos.clear()
    circos.par("gap.degree"=1, start.degree = 90)
    circos.initializeWithIdeogram(cytoband.df)
    circos.genomicLink(bedL, bedR, col=ipet.col)
    pushViewport(viewport(x=0.8, y=0.05, width=0.4, height=0.05, name="Leg", just=c("center","bottom")))
    grid.draw(lgd_list_horizontal);
    upViewport(0);
    title(paste0(cell,"Chia-PET_chiadrop_TSS_ecSE1_","_Norm",group,": ",nrow(itxs)))
    dev.off()
}

#links for circos within ecDNA (ChIA-PET data)
itx = read.table("/net/nwgc/vol1/nobackup/nocleanup/tungch/test/diffloop/trans_colo/ec_seperated_connected/normalized_byMean/SE_annotation_with_ec/PC3DM+.drug_vs_mock.diffDAitxTrans.chiapet_chr_basedOnChiadrop_gene.with_ec_anchors.ecSEannotated.tsv", stringsAsFactors=F, header = TRUE)

# mock
itx2 = itx[order(itx$norm_mock),]
table(itx2$SErank_L == 1 )
which(itx2$SErank_L == 1 )
mock = subset(itx2, norm_mock != 0 & SErank_L == 1, select = c(2:7,12))
dim(mock)

# HD
itx2 = itx[order(itx$norm_drug),]
HD = subset(itx2, norm_drug != 0 & SErank_L == 1, select = c(2:7,13))
dim(HD)

process_samples(mock,"PC3","mock",2.5,3,"#f40000")
process_samples(HD,"PC3","HD",2.5,3,"#f40000")
