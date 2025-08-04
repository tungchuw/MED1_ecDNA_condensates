##################################
# GSEA analysis 
# Updated: 12.19.2023
# Author: Chia-Hao T.
###################################

###################################

library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(ggpubr)
library(aplot)

###################################

cellname="COLO320DM"

## Generate plot theme.
plot_theme <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 10),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))


# Load mastertable
rnaseq_se <- read.delim("./data/misc/COLO320DM/rnaseq/COLO320DM.AT4-48.rnaSeq.DESeq2.MasterTable.annotated.tsv") 
rnaseq_se$ENS_id <- gsub("\\..*","", rnaseq_se$ENS)
rnaseq_se_unique <- rnaseq_se[!duplicated(rnaseq_se$ENS_id),]

# GSEA analysis of Combo_SE1to4 prerank
sortdf_Combo_SE1to4<-rnaseq_se_unique[order(rnaseq_se_unique$log2FoldChange_Combo_SE1to4, decreasing = T),]
head(sortdf_Combo_SE1to4)
prerank_Combo_SE1to4 = sortdf_Combo_SE1to4$log2FoldChange_Combo_SE1to4
head(prerank_Combo_SE1to4)
names(prerank_Combo_SE1to4) <- sortdf_Combo_SE1to4$ENS_id
head(prerank_Combo_SE1to4)

gsea_BP_Combo_SE1to4 <- gseGO(
  prerank_Combo_SE1to4,
  ont = "BP",   
  OrgDb = "org.Hs.eg.db", 
  keyType = "ENSEMBL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
)

# Save all ontolgies sorted by adjusted pvalues
write.table(gsea_BP_Combo_SE1to4, file=paste0("Summary_GSEA_analysis_COLO320DM_BP_Combo_SE1to4.tsv"), sep='\t', quote=F, col.names=T, row.names=F)

#plot gsea
a <- gseaplot2(gsea_BP_Combo_SE1to4, geneSetID = which(gsea_BP_Combo_SE1to4$Description == "DNA replication"), title = paste0(cellname,"_comboSE_","DNA replication","_Padj = ",gsea_BP_Combo_SE1to4$p.adjust[which(gsea_BP_Combo_SE1to4$Description == "DNA replication")],"_n = ",gsea_BP_Combo_SE1to4$setSize[which(gsea_BP_Combo_SE1to4$Description == "DNA replication")]))

pdf(file=paste0(cellname,'.gsea_plots.pdf'), width=8, height=6)
print(a)
dev.off()


