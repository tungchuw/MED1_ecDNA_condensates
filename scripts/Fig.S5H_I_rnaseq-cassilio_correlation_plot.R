
library(pheatmap)
library(dplyr)

rna_PC3 = read.delim(paste0('/net/nwgc/vol1/sharing/Wei_Lab/proj-SuperEnhancer/rnaseq/transfection/deSeq/AT4-31.rnaSeq.DE_emptyVec.MasterTable.tsv'), stringsAsFactors=F)
rna_COLO = read.delim(paste0('/net/nwgc/vol1/nobackup/nocleanup/tungch/RNA_AT4-48/deSeq/without_ctrl/RNA_AT4-48.rnaSeq.DE_emptyVec_woCtrlSE.MasterTable.tsv'), stringsAsFactors=F)

normCounts_rna_PC3 = rna_PC3[,c(2:19)]
normCounts_rna_COLO = rna_COLO[,c(2:19)]

cormat_normCounts_rna_PC3 = cor(normCounts_rna_PC3) 
cormat_normCounts_rna_COLO = cor(normCounts_rna_COLO) 

# with numbers
pdf(paste0('PC3_casilio_RNAseq_correlation_plot.pdf'), width=6, height=6)
pheatmap(cormat_normCounts_rna_PC3, 
        display_numbers=T, 
        border_color = '#c0c0c0c0', 
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        # number_color = "black",
        fontsize_number=6, 
        fontsize=8, 
        main=paste0('PC3_casilio_RNAseq_correlation'))
dev.off()


# with numbers
pdf(paste0('COLO320DM_casilio_RNAseq_correlation_plot.pdf'), width=6, height=6)
pheatmap(cormat_normCounts_rna_COLO, 
        display_numbers=T, 
        border_color = '#c0c0c0c0', 
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        # number_color = "black",
        fontsize_number=6, 
        fontsize=8, 
        main=paste0('COLO320DM_casilio_RNAseq_correlation'))
dev.off()


