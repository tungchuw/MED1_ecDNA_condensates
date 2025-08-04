load('/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/1x51_chipseq_processed/diffbind/HD2/B171/diffbind3.8.4.B171.RData')

library(DiffBind)
library(dplyr)
library("ggplot2")
library(pheatmap)

# with numbers
pdf(paste0('B171_HD_MED1ChIPseq_correlation_plot.pdf'), width=6, height=6)
# mycolor = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
# ph = dba.plotHeatmap(db.obj)
ph_1 = ph[c("A0088","A0089","A0090","A0091","A0092","A0093"),c("A0088","A0089","A0090","A0091","A0092","A0093")]
pheatmap(ph_1, 
        display_numbers=T, 
        border_color = '#c0c0c0c0', 
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        # number_color = "black",
        fontsize_number=12, 
        fontsize=13, 
        main=paste0('B171_HD_MED1ChIPseq_correlation'))
dev.off()
