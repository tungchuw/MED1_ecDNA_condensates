load('/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/1x51_chipseq_processed/diffbind/HD2/COLO320DM/diffbind3.8.4.COLO320DM.RData')

library(DiffBind)
library(dplyr)
library("ggplot2")
library(pheatmap)

# with numbers
pdf(paste0('COLO320DM_HD_MED1ChIPseq_correlation_plot.pdf'), width=6, height=6)
# mycolor = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
ph = dba.plotHeatmap(db.obj)
ph_1 = ph[c("A0076","A0081","A0077","A0082","A0078","A0083","A0079","A0084","A0080","A0086"),c("A0076","A0081","A0077","A0082","A0078","A0083","A0079","A0084","A0080","A0086")]
pheatmap(ph_1, 
        display_numbers=T, 
        border_color = '#c0c0c0c0', 
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        # number_color = "black",
        fontsize_number=8, 
        fontsize=13, 
        main=paste0('COLO320DM_HD_MED1ChIPseq_correlation'))
dev.off()
