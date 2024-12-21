if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
library(dplyr)
library(readr)
library(data.table)
library(pheatmap)
library(tidyr)

# Have calculated the correlation matrix in Plink2 and have these saved to the LD directory on UKB RAP
# dx download /Output/regenie/regenie_LD/
# Read in the priority variants 
# dx download /Output/gene_VCF_variants/variants/pri_variants/chrpos/chrposrefalt_gene_canonical/FADS1_priority_annot_chrposrefalt.txt
# dx download /Input/FADS_cluster_UKB_pVCF.tsv

# read in FADS file 
fads <- read.table("FADS_cluster_UKB_pVCF.tsv", sep = "\t", header =T)

plot_heatmaps <- function(gene) {
    print(gene)
    print(paste0("Reading in correlation matrix calculated by REGENIE: mdd_", gene, "_step2_BT_LD_cor.corr"))
# Read in the correlation matrix
cor <- fread(paste0("mdd_", gene, "_step2_BT_LD_cor.corr")) %>% as.data.frame()
# Read in the variant names 
cor_names <- readr::read_lines(paste0("mdd_", gene, "_step2_BT_LD_cor.corr.snplist"))
colnames(cor) <- cor_names
row.names(cor) <- cor_names

# Priority variants 
priority <- readr::read_lines(paste0(gene, "_priority_annot_chrposrefalt.txt"))
priority <- gsub("chr", "", priority)
priority <- gsub("_", ":", priority)

# Filter the matrix to the priority variants
filtered_cor <- cor[row.names(cor) %in% priority, colnames(cor) %in% priority]
png(paste0(gene, "_priority_heatmap.png"), width = 8, height = 6, units = "in", res = 300)
pheatmap(filtered_cor, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

zoomed_cor <- filtered_cor
zoomed_cor$Row <- rownames(zoomed_cor)
cor_long <- zoomed_cor %>%
  pivot_longer(cols = -Row, names_to = "Var2", values_to = "Correlation") %>%
  rename(Var1 = Row) 
variants_highercor <- cor_long %>% filter(Correlation > 0.25 & Var1 != Var2) %>% pull(Var1)
if (length(variants_highercor > 0)) {
    print("Plotting the variants with more than 0.25 correlation")
    # Filter to those with a correlation > 0.25
higher_cor_mat <- filtered_cor[row.names(filtered_cor) %in% variants_highercor, colnames(filtered_cor) %in% variants_highercor]
png(paste0(gene, "_higher_corr_heatmap.png"), width = 8, height = 6, units = "in", res = 300)
pheatmap(higher_cor_mat)
dev.off()
} else {
    print("No variants with > 0.25 correlation")
}
}

for(i in 1:nrow(fads)) {
    gene_name <- fads$hgnc_symbol[i]
    print(gene_name)
    plot_heatmaps(gene_name)
}

