setwd("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Intermediate_data_files/vep_res/AC_overzero_2712/")
library(gridExtra)
library(grid)
library(ggplot2)

FADS <- read.table("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Intermediate_data_files/FADS_cluster_UKB_pVCF.tsv", sep = "\t", header = T)
for (gene in FADS$hgnc_symbol) {
  print(gene)
  
# FADS1 
imgA=png::readPNG(paste0(gene, "_priority_plots/UpSet_", gene, ".png"))
imgB=png::readPNG(paste0(gene, "_priority_plots/CSQ_priority", gene, ".png"))
imgC=png::readPNG(paste0(gene, "_priority_plots/AC_hist_prioritised", gene, ".png"))

# Convert images to grobs
grobA <- rasterGrob(imgA, interpolate = TRUE)
grobB <- rasterGrob(imgB, interpolate = TRUE)
grobC <- rasterGrob(imgC, interpolate = TRUE)

add_label <- function(grob, label) {
  label_grob <- textGrob(
    label, x = unit(0, "npc") + unit(5, "mm"), y = unit(1, "npc") - unit(5, "mm"),
    hjust = 0, vjust = 1, gp = gpar(fontsize = 12, fontface = "bold")
  )
  # Overlay label on the plot
  gTree(children = gList(grob, label_grob))
}

# Add labels to each plot
labeledA <- add_label(grobA, "A)")
labeledB <- add_label(grobB, "B)")
labeledC <- add_label(grobC, "C)")

# Create a title
title <- textGrob(
  paste0(gene, ": Prioritised variants"), 
  gp = gpar(fontsize = 17, fontface = "bold"), 
  x = 0.5, hjust = 0.5
)

plot_dir <- "/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Plots/"
# Arrange title and plots
pdf(paste0(plot_dir, gene, "_vep_priority_arranged.pdf"), width = 8.27, height = 11.69) # A4 dimensions
grid.arrange(
  title, labeledA, labeledB, labeledC,
  ncol = 1, heights = c(0.04, 0.32, 0.32, 0.32) # Adjust heights for title and plots
)
dev.off()
}
