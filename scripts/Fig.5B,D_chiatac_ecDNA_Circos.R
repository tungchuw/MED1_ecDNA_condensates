# scripts/Fig.5B_ecDNA_Circos.R
# Generate ecDNA-only circos plots with normalized iPET signal tracks

suppressMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(gtools)
  library(zoo)
  library(gridBase)
  library(dplyr)
  library(stringr)
  library(gplots)
})

# Load cytoband info for hg38
inclCHR <- paste0("chr", c(1:22, "X"))
cytoband.df <- subset(
  readRDS(system.file("extdata", "cytoband_list.rds", package = "circlize"))[["hg38"]],
  V1 %in% inclCHR
)

# Read interaction table
itx <- read.table(
  "./data/misc/PC3DM/chiatac/PC3DM+.trans-loops_normalized_annotated.tsv",
  header = TRUE, stringsAsFactors = FALSE
)

# Filter and extract vector, SE1, combo samples with SErank_L == 1
# vector
vector = subset(itx, norm_AT4.40.19m20.vector_iPET != 0 & SErank_L == 1, select = c(2:7, 11))
vector = vector[order(vector$norm_AT4.40.19m20.vector_iPET), ]
# SE1
SE1 = subset(itx, norm_AT4.40.1m2.SE1_iPET != 0 & SErank_L == 1, select = c(2:7, 12))
SE1 = SE1[order(SE1$norm_AT4.40.1m2.SE1_iPET), ]
# combo
combo = subset(itx, norm_AT4.40.9m10.combo_SE_iPET != 0 & SErank_L == 1, select = c(2:7, 13))
combo = combo[order(combo$norm_AT4.40.9m10.combo_SE_iPET), ]

# Plot parameters
min_val <- 2.5
max_val <- 6

# Main plot function
generate_circos <- function(itxs, cell, group, min_val, max_val, max_color) {
  message(paste0("Plotting ", group, " (n = ", nrow(itxs), ")"))
  
  bedL <- itxs[, 1:3]
  bedR <- itxs[, 4:6]
  ipet <- itxs[, 7]

  # Color ramp setup
  breaks <- c(0, 1:10, max(ipet) + 1)
  interval <- cut(ipet, breaks = breaks, right = FALSE)
  print(table(interval))

  col_fun <- colorRamp2(c(min_val, (min_val + max_val)/2, max_val), 
                        c("#E7E7E7", "orange", max_color), transparency = 0.5)
  ipet.col <- col_fun(ipet)

  # Legend
  lgd_link <- Legend(
    at = c(min_val, (min_val + max_val)/2, max_val),
    col_fun = col_fun,
    title_position = "topcenter",
    title = "Normalized iPET count",
    direction = "horizontal"
  )
  
  # Begin plotting
  pdf(paste0(cell, "_ChIATAC_ecSE1_", group, "_min_", min_val, "_max_", max_val, "_scale_color_", max_color, "_updated.pdf"), 
      width = 8, height = 8, pointsize = 15)
  circos.clear()
  circos.par("gap.degree" = 1, start.degree = 90)
  circos.initializeWithIdeogram(cytoband.df)
  circos.genomicLink(bedL, bedR, col = ipet.col)

  pushViewport(viewport(x = 0.8, y = 0.05, width = 0.4, height = 0.05, just = c("center", "bottom")))
  grid.draw(packLegend(lgd_link))
  upViewport(0)
  dev.off()
}

# Run for all groups
generate_circos(vector, "PC3", "vector", min_val, max_val, "#f40000")
generate_circos(SE1, "PC3", "SE1", min_val, max_val, "#f40000")
generate_circos(combo, "PC3", "combo", min_val, max_val, "#f40000")

message("Finished generating ecDNA circos plots.")
