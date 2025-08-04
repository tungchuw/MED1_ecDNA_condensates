load('diffbind3.8.4_pc3dm+.span550.RData')

library(DiffBind)
library(dplyr)
library("ggplot2")
library(pheatmap)

# with numbers
pdf(paste0('PC3_HD_MED1ChIPseq_correlation_plot.pdf'), width=6, height=6)
# mycolor = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
# ph = dba.plotHeatmap(db.obj)
ph_1 = ph[c("A0069","A0072","A0070","A0073","A0071","A0074"),c("A0069","A0072","A0070","A0073","A0071","A0074")]
pheatmap(ph_1, 
        display_numbers=T, 
        border_color = '#c0c0c0c0', 
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        # number_color = "black",
        fontsize_number=12, 
        fontsize=13, 
        main=paste0('PC3_HD_MED1ChIPseq_correlation'))
dev.off()

colnames(ph)

