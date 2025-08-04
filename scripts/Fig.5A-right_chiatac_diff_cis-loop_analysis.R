# scripts/Fig.5A_right_diff_cis-iPET_analysis.R
# Main pipeline: iPET differential interaction analysis, annotation, and export

suppressMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(pheatmap)
  library(ggrepel)
  library(ggbreak)
  library(scales)
  library(ggbeeswarm)
})

# === CONFIGURATION ===
cellname <- "PC3DM+"
itxsuf <- ".DA.itx"

# Directory paths
base_dir <- "./data/misc/"
cell_dir <- "./data/misc/PC3DM/"
script_dir <- "./scripts/"

# External tools
bedtools <- paste0("singularity run ./singularity/chipseqtools.sif bedtools intersect -wo -a ")
bedmerge <- paste0("singularity run ./singularity/chipseqtools.sif bedtools merge -i ")

# Input files
sample_file <- file.path(cell_dir, paste0("chiatac/",cellname, "_all_itxlist.txt"))
samplePhenotype <- read.table(sample_file)
colnames(samplePhenotype) <- c("name","condition", "rep")
itxlibs <- samplePhenotype$name

# Regions
ecDNAfile <- file.path(cell_dir, paste0(cellname, ".ecDNA.bed"))
peakfile <- file.path(cell_dir, paste0(cellname, ".Pol2.narrowPeak"))
roseBed <- file.path(cell_dir, paste0(cellname , ".SErose.bed"))
med1peak <- file.path(cell_dir, paste0(cellname, ".MED1overlap_peak.bed"))
genebed <- file.path(base_dir, "hg38.Genes.itxAnnotation.bed")
perlannot <- paste0("/usr/bin/perl ", file.path(script_dir, "itxAnnotation.4peaks.pl"))

# Output directories
out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE)

# Step 1: Load interaction tables
itxs <- list()
headerset <- c('ChrL', 'StartL', 'EndL', 'ChrR', 'StartR', 'EndR', 'iPET')
for (i in itxlibs) {
  f1 <- file.path(cell_dir, "chiatac/", paste0(i, itxsuf))
  message("Reading: ", f1)
  m1 <- fread(f1, header = FALSE)[, 1:7]
  colnames(m1) <- headerset
  itxs[[i]] <- m1
}

# Step 2: Create anchor nodes by merging all ends of iPETs
allitx = bind_rows(itxs) #combine all itx
anchors1 = allitx[,1:3]
anchors2 = allitx[,4:6]
colnames(anchors1) = colnames(anchors2) = c('CHR','Start','End')
anchors = rbind(anchors1, anchors2)
rm(anchors1,anchors2)

# Save and merge anchors
tmpbed <- paste0(cellname, ".tmp.bed")
tmpout <- paste0(cellname, ".tmp.out")
write.table(anchors, tmpbed, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
system(paste0("sort -k1,1 -k2,2n ", tmpbed, " | ", bedmerge, " stdin > ", tmpout))
anchors <- fread(tmpout)
colnames(anchors) <- c("CHR", "Start", "End")
anchors$ID <- paste0(anchors$CHR, ":", anchors$Start, "-", anchors$End)
anchorfile <- paste0(cellname, ".anchors.bed")
write.table(anchors, anchorfile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Annotate ecDNA status
system(paste0(bedtools, anchorfile, " -b ", ecDNAfile, " > ", tmpout))
overlap <- read.table(tmpout)
atype <- rep("ch", nrow(anchors))
atype[anchors$ID %in% overlap$V4] <- "ec"
anchors$locus <- atype

# Step 3: Define function to map interactions to anchors
getAnchorIndex <- function(itx1, libID) {
  n <- nrow(itx1)
  tmpbed <- "tmp.bed"; tmpout <- "tmp.txt"

  anchorL <- numeric(n)
  anchorR <- numeric(n)

  write.table(cbind(itx1[, 1:3], 1:n), tmpbed, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  system(paste0(bedtools, tmpbed, " -b ", anchorfile, " > ", tmpout))
  anchorL[fread(tmpout)$V4] <- fread(tmpout)$V8

  write.table(cbind(itx1[, 4:6], 1:n), tmpbed, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  system(paste0(bedtools, tmpbed, " -b ", anchorfile, " > ", tmpout))
  anchorR[fread(tmpout)$V4] <- fread(tmpout)$V8

  anchorLab <- paste0(pmin(anchorL, anchorR), "_", pmax(anchorL, anchorR))
  idPETs <- data.frame(pairID = anchorLab, iPET = itx1[, 7], stringsAsFactors = FALSE)
  idPETs <- aggregate(iPET ~ pairID, data = idPETs, FUN = sum)
  colnames(idPETs)[2] <- paste0(libID, "_iPET")
  return(idPETs)
}

# Step 4: Generate raw counts per anchor pair
anchor.list <- lapply(names(itxs), function(i) getAnchorIndex(itxs[[i]], i))
names(anchor.list) <- names(itxs)
rawcounts <- Reduce(function(x, y) merge(x, y, by = "pairID", all = TRUE), anchor.list)
rawcounts[is.na(rawcounts)] <- 0
rownames(rawcounts) <- rawcounts$pairID
rawcounts <- rawcounts[, -1]

# Step 5: Normalize and annotate
normFactor <- colMeans(rawcounts)
normMat <- sweep(rawcounts, 2, normFactor, FUN = "/")
colnames(normMat) <- paste0("norm_", colnames(normMat))
DEtable <- cbind(rawcounts, normMat)
DEtable$itxID <- rownames(DEtable)

# Step 6: Generate coordinate and annotation
nodeLR <- strsplit(DEtable$itxID, "_")
get_coords <- function(vec) {
  ch <- sapply(strsplit(vec, ":"), "[", 1)
  pos <- sapply(strsplit(vec, ":"), "[", 2)
  pos_split <- strsplit(pos, "-")
  list(chr = ch, start = sapply(pos_split, "[", 1), end = sapply(pos_split, "[", 2))
}
coordL <- get_coords(sapply(nodeLR, "[", 1))
coordR <- get_coords(sapply(nodeLR, "[", 2))
tmpTable <- data.frame(chrL = coordL$chr, startL = coordL$start, endL = coordL$end,
                       chrR = coordR$chr, startR = coordR$start, endR = coordR$end,
                       p = 0, ipet = 1, qval = 0)
tmpname <- paste0(cellname, ".tmp.txt")
write.table(tmpTable, tmpname, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
system(paste(perlannot, tmpname, genebed, peakfile, ecDNAfile, roseBed, med1peak, cellname))

# Step 7: Load annotation results
itx <- fread(paste0(cellname, ".annotated_itx.txt"), header = TRUE)
uproL <- sapply(strsplit(itx$TSS_L, ";"), function(x) paste(unique(x), collapse = ";"))
uproR <- sapply(strsplit(itx$TSS_R, ";"), function(x) paste(unique(x), collapse = ";"))

# Step 8: Compute SErank
get_SErank <- function(df, anchor_cols) {
  write.table(cbind(df[, anchor_cols], 1:nrow(df)), tmpname, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  system(paste0(bedtools, tmpname, " -b ", roseBed, " > ", tmpout))
  o <- read.table(tmpout)
  ranks <- numeric(nrow(df))
  ranks[o$V4] <- o$V8
  return(ranks)
}

SErankL <- get_SErank(tmpTable, 1:3)
SErankR <- get_SErank(tmpTable, 4:6)

# Step 9: Final master table
mastertab <- data.frame(DEtable, ecDNA = itx$ecDNA, PGI = itx$FinalCode,
                        SErank_L = SErankL, SErank_R = SErankR,
                        MED1 = itx$MED1, TSS_L = uproL, TSS_R = uproR,
                        chr_L = itx$chrL, start_L = itx$startL, end_L = itx$endL,
                        chr_R = itx$chrR, start_R = itx$startR, end_R = itx$endR)

# Save output
out_file <- file.path(out_dir, paste0(cellname, ".node_interaction_master.tsv"))
write.table(mastertab, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

message("iPET interaction analysis completed and saved to:", out_file)


# Step 10: Summary stats and violin plot for LR ecDNA interactions

# Count of LR ecDNA iPETs across libraries
vec_col <- grep("norm.*vector_iPET", colnames(mastertab), value=TRUE)
se1_col <- grep("norm.*SE1_iPET", colnames(mastertab), value=TRUE)
combo_col <- grep("norm.*combo_SE_iPET", colnames(mastertab), value=TRUE)

LR_mask <- mastertab$ecDNA == "LR"
nonzero_vec <- mastertab[[vec_col]] != 0 & LR_mask
nonzero_se1 <- mastertab[[se1_col]] != 0 & LR_mask
nonzero_combo <- mastertab[[combo_col]] != 0 & LR_mask

Values_vec <- mastertab[[vec_col]][nonzero_vec]
Values_se1 <- mastertab[[se1_col]][nonzero_se1]
Values_combo <- mastertab[[combo_col]][nonzero_combo]

# Wilcoxon tests
pval_vec_vs_combo <- wilcox.test(Values_combo, Values_vec, alternative = "less")$p.value
pval_vec_vs_se1 <- wilcox.test(Values_se1, Values_vec, alternative = "less")$p.value

# Violin plot prep
group_labels <- c(rep("vector", length(Values_vec)), 
                  rep("SE1", length(Values_se1)), 
                  rep("combo_SE", length(Values_combo)))

cp.df <- data.frame(normFreq = c(Values_vec, Values_se1, Values_combo), 
                    cond = factor(group_labels, levels = c("vector", "SE1", "combo_SE")))

# Summary stats (optional printout)
summary(Values_vec)
summary(Values_se1)
summary(Values_combo)

# Violin plot output
pdf(paste0(cellname, 'MergedFastq.compare_cis_non-pairwise_SE1.pdf'), width = 6, height = 6)
ggplot(cp.df, aes(x = cond, y = normFreq, fill = cond)) + 
  geom_violin(width = 1, bw = 0.05) +
  scale_x_discrete(labels = c(paste0("vector (n = ", length(Values_vec), ")"), 
                              paste0("SE1 (n = ", length(Values_se1), ")"), 
                              paste0("combo_SE (n = ", length(Values_combo), ")"))) +
  xlab("") +
  scale_fill_manual(values = c("vector" = "#7fbf7b", "SE1" = "#af8dc3", "combo_SE" = "#af8dc3")) +
  geom_boxplot(width = 0.03, color = "black", alpha = 10) +
  scale_y_continuous(trans = 'log10', breaks = c(0, 1, 10, 100), guide = "axis_logticks") +
  theme(legend.position = "bottom", panel.spacing.x = unit(1.5, "lines")) +
  theme_classic()
dev.off()

message("Statistical testing and violin plot generated.")

# clean

# system("rm tmp* PC3*")
