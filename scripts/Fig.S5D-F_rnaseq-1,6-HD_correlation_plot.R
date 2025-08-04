
library(pheatmap)
library(dplyr)

rna_PC3 = read.delim(paste0('/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/rnaseq/HDtreatment/PC3DM+/withChiadrop/PC3DM+.RNAseq_HD2.MasterTable.annotated.tsv'), stringsAsFactors=F)
rna_COLO = read.delim(paste0('/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/rnaseq/HDtreatment/COLO320DM/withChiadrop/COLO320DM.RNAseq_HD2.MasterTable.annotated.tsv'), stringsAsFactors=F)
rna_B171 = read.delim(paste0('/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/rnaseq/HDtreatment/B171/withChiadrop/B171.RNAseq_HD2.MasterTable.annotated.tsv'), stringsAsFactors=F)

normCounts_rna_PC3 = subset(rna_PC3,select = c("ARP0038","ARP0041","ARP0039","ARP0042","ARP0040","ARP0043"))
normCounts_rna_COLO = subset(rna_COLO,select = c("ARP0044","ARP0048","ARP0045","ARP0049","ARP0047","ARP0051"))
normCounts_rna_B171 = subset(rna_B171,select = c("ARP0052","ARP0053","ARP0054","ARP0055","ARP0056","ARP0057"))

cormat_normCounts_rna_PC3 = cor(normCounts_rna_PC3) 
cormat_normCounts_rna_COLO = cor(normCounts_rna_COLO) 
cormat_normCounts_rna_B171 = cor(normCounts_rna_B171) 

# with numbers
pdf(paste0('PC3_HD_RNAseq_correlation_plot.pdf'), width=6, height=6)
pheatmap(cormat_normCounts_rna_PC3, 
        display_numbers=T, 
        border_color = '#c0c0c0c0', 
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        # number_color = "black",
        fontsize_number=8, 
        fontsize=13, 
        main=paste0('PC3_HD_RNAseq_correlation'))
dev.off()

# with numbers
pdf(paste0('COLO320DM_HD_RNAseq_correlation_plot_remove2h.pdf'), width=6, height=6)
pheatmap(cormat_normCounts_rna_COLO, 
        display_numbers=T, 
        border_color = '#c0c0c0c0', 
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        # number_color = "black",
        fontsize_number=8, 
        fontsize=13, 
        main=paste0('COLO320DM_HD_RNAseq_correlation'))
dev.off()

# with numbers
pdf(paste0('B171_HD_RNAseq_correlation_plot.pdf'), width=6, height=6)
pheatmap(cormat_normCounts_rna_B171, 
        display_numbers=T, 
        border_color = '#c0c0c0c0', 
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        # number_color = "black",
        fontsize_number=8, 
        fontsize=13, 
        main=paste0('B171_HD_RNAseq_correlation'))
dev.off()
