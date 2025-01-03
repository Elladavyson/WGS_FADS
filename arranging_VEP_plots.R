
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

# Lolliplot subplots for different masks in each gene

masks <- c("Mask1.0.01", "Mask1.singleton", "Mask2.0.01", "Mask2.singleton",
           "Mask3.0.01", "Mask3.singleton", "Mask4.0.01", "Mask4.singleton",
           "Mask5.0.01", "Mask5.singleton")

for (gene in FADS$hgnc_symbol) {
  print(gene)
  # FADS1 
  imgA=png::readPNG(paste0("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Plots/", gene, "_lolliplot_variants_all.png"))
  grobA <- rasterGrob(imgA)
  assign(paste0(gene, "_lolliplots"), imgA)
  assign(paste0(gene, "_lolliplots_grob"), grobA)
}

for (gene in FADS$hgnc_symbol) {
  print(gene)
  mask_grobs <- list()
for (mask in masks[grepl("0.01", masks)]) {
  png_path <- paste0("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Plots/", gene, "_", mask, "_lolliplot_variants_pri.png")
if (file.exists(png_path)) {
  masklolliplot <- png::readPNG(paste0("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Plots/", gene, "_", mask, "_lolliplot_variants_pri.png"))
  mask_grobs[[mask]] <- rasterGrob(masklolliplot)
} else {
  print(paste0("File not found:", png_path))
}
}
  if (length(mask_grobs) > 0) {
    labeled_grobs <- mapply(add_label, mask_grobs, LETTERS[seq_along(mask_grobs)], SIMPLIFY = FALSE)
    pdf_path <- paste0(plot_dir, gene, "_lolliplots_all0.01masks.pdf")
    pdf(pdf_path, width = 8.27, height = 11.69)
    grid.arrange(grobs = labeled_grobs, ncol = 1, heights = rep(0.2, length(labeled_grobs)))
    dev.off()
    print(paste0("Saved PDF for gene 0.01 masks: ", gene, " at ", pdf_path))
  }  else {
      print(paste0("No valid plots found for gene: ", gene))
  }
  mask_grobs <- list()
  for (mask in masks[grepl("singleton", masks)]) {
    png_path <- paste0("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Plots/", gene, "_", mask, "_lolliplot_variants_pri.png")
    if (file.exists(png_path)) {
      masklolliplot <- png::readPNG(paste0("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Plots/", gene, "_", mask, "_lolliplot_variants_pri.png"))
      mask_grobs[[mask]] <- rasterGrob(masklolliplot)
    } else {
      print(paste0("File not found:", png_path))
    }
  }
  if (length(mask_grobs) > 0) {
    labeled_grobs <- mapply(add_label, mask_grobs, LETTERS[seq_along(mask_grobs)], SIMPLIFY = FALSE)
    pdf_path <- paste0(plot_dir, gene, "_lolliplots_allsingletonmasks.pdf")
    pdf(pdf_path, width = 8.27, height = 11.69)
    grid.arrange(grobs = labeled_grobs, ncol = 1, heights = rep(0.2, length(labeled_grobs)))
    dev.off()
    print(paste0("Saved PDF for gene singleton masks: ", gene, " at ", pdf_path))
  }  else {
    print(paste0("No valid plots found for gene: ", gene))
  }
  }

# All prioritised variants FADS1-3
  plot_dir <- "/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Plots/"
  # Arrange title and plots
  pdf(paste0(plot_dir, "FADS1_3_lolliplots_all.pdf"), width = 8.27, height = 11.69) # A4 dimensions
  grid.arrange(add_label(FADS1_lolliplots_grob, "A)"),add_label(FADS2_lolliplots_grob, "B)"), add_label(FADS3_lolliplots_grob, "C)"),
    ncol = 1, heights = c(0.34, 0.33, 0.33) # Adjust heights for title and plots
  )
  dev.off()

# All prioritised variants FEN1, MYRF, TMEM258
  pdf(paste0(plot_dir, "FEN1_MYRF_TMEM258_lolliplots_all.pdf"), width = 8.27, height = 11.69) # A4 dimensions
  grid.arrange(add_label(FEN1_lolliplots_grob, "A)"),add_label(MYRF_lolliplots_grob, "B)"), add_label(TMEM258_lolliplots_grob, "C)"),
               ncol = 1, heights = c(0.34, 0.33, 0.33) # Adjust heights for title and plots
  )
dev.off()
  
# ALL prioritised variants for all genes

pdf(paste0(plot_dir, "all_lolliplots.pdf"), width = 8.27, height = 11.69) # A4 dimensions
grid.arrange(add_label(FADS1_lolliplots_grob, "A)"),
             add_label(FADS2_lolliplots_grob, "B)"),
             add_label(FADS3_lolliplots_grob, "C)"),
             add_label(FEN1_lolliplots_grob, "D)"),
             add_label(MYRF_lolliplots_grob, "E)"), 
             add_label(TMEM258_lolliplots_grob, "F)"),
             ncol = 1, heights = c(0.16, 0.17, 0.16, 0.17, 0.17, 0.17) # Adjust heights for title and plots
)
dev.off()
  

# Plotting the beta distribution for weights in R 

alpha = 1
beta = 25
x <- seq(0,1, length.out = 100)
y <- dbeta(x, shape = alpha, shape2= beta)
beta_df <- data.frame(x = x, y = y)
beta_dist <- ggplot(data = beta_df, aes(x=x, y=y)) + 
  geom_line(color = "blue", size = 1) +
  labs(x = "Minor Allele Frequency", y = "Weight", title=
  expression("Beta Distribution " ~ Beta(MAF[i], 1, 25)))+ theme_bw()

ggsave(paste0(plot_dir, "beta_weights_distribution.png"), plot = beta_dist, width = 6, height = 6, device='png', dpi=300)

## Variant summaries per gene 
for (gene in FADS$hgnc_symbol) {
  print(gene)
  carrier_summary=png::readPNG(paste0("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Plots/", gene, "_pri_carriers_summary.png"))
  variant_summary=png::readPNG(paste0("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Plots/", gene, "_pri_variant_summary.png"))
  grobcarriers <- rasterGrob(carrier_summary)
  grobvariants <- rasterGrob(variant_summary)
  assign(paste0(gene, "_carriersummary_grob"), grobcarriers)
  assign(paste0(gene, "_variantsummary_grob"), grobvariants)
}

pdf(paste0(plot_dir, "all_carrier_summaries.pdf"), width = 8.27, height = 11.69) # A4 dimensions
grid.arrange(add_label(FADS1_carriersummary_grob, "A) FADS1"),
             add_label(FADS2_carriersummary_grob, "B) FADS2"),
             add_label(FADS3_carriersummary_grob, "C) FADS3"),
             add_label(FEN1_carriersummary_grob, "D) FEN1"),
             add_label(MYRF_carriersummary_grob, "E) MYRF"), 
             add_label(TMEM258_carriersummary_grob, "F) TMEM258"),
             ncol = 2, nrow = 3, heights = c(0.3, 0.3, 0.3) # Adjust heights for title and plots
)
dev.off()

pdf(paste0(plot_dir, "all_variant_summaries.pdf"), width = 8.27, height = 11.69) # A4 dimensions
grid.arrange(add_label(FADS1_variantsummary_grob, "A) FADS1"),
             add_label(FADS2_variantsummary_grob, "B) FADS2"),
             add_label(FADS3_variantsummary_grob, "C) FADS3"),
             add_label(FEN1_variantsummary_grob, "D) FEN1"),
             add_label(MYRF_variantsummary_grob, "E) MYRF"), 
             add_label(TMEM258_variantsummary_grob, "F) TMEM258"),
             ncol = 2, nrow = 3, heights = c(0.3, 0.3, 0.3) # Adjust heights for title and plots
)
dev.off()



