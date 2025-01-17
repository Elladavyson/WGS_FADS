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
# dx download /Output/priority_LD
# dx download /Output/genotypes/genotype_summary/all_priority/*_carriers_genotypes.tsv
# dx download /Input/FADS_cluster_UKB_pVCF.tsv

# read in FADS file 
fads <- read.table("FADS_cluster_UKB_pVCF.tsv", sep = "\t", header =T)

plot_heatmaps <- function(gene) {
    print(gene)
    print(paste0("Reading in correlation matrix calculated by PLINK:", gene, "_priority_LD.ld"))
# Read in the LD file from PLINK 
ld <- read.table(paste0(gene, "_priority_LD.ld"), header =T) %>% as.data.frame()
print("Reading in carriers genotypes to check the higher LD variants")
genotypes <- read.table(paste0(gene, "_carriers_genotypes.tsv"), sep ="\t", header =T)
unique_variants <- unique(c(ld$SNP_A, ld$SNP_B))
# LD matrix of all priority variants 
ld_matrix <- ld %>%
  select(SNP_A, SNP_B, R2) %>%
  pivot_wider(names_from = SNP_B, values_from = R2, values_fill = 0) %>% as.data.frame()
row.names(ld_matrix) <- ld_matrix$SNP_A
ld_matrix <- ld_matrix %>% select(-SNP_A)
ld_matrix <- ld_matrix[unique_variants, unique_variants, drop = FALSE]

# Priority variants 
png(paste0(gene, "_priority_heatmap.png"), width = 8, height = 6, units = "in", res = 300)
pheatmap(as.matrix(ld_matrix), 
cluster_cols = FALSE,
cluster_rows = FALSE,
 show_rownames = FALSE,
 show_colnames= FALSE,
 display_numbers = FALSE,
 breaks = seq(0,1, length.out = 101),
 main = paste0(gene, ": Prioritised variants R2 (n =", length(unique(ld$SNP_A)), ")"))
dev.off()

# Proportion of distinct variant pairs variants under a certain correlation 
ld_distinct_pairs <- ld %>% filter(SNP_A != SNP_B)
ld_over0.01 <- ld_distinct_pairs %>% filter(R2 > 0.01)
max(ld_distinct_pairs$R2) # The maximum R2 for the prioritised variants in this gene 
if (ld_distinct_pairs %>% filter(R2==1) %>% nrow()> 0) {
ld_1 <- ld_distinct_pairs %>% filter(R2 == 1) %>% pull(SNP_A) %>% unique() 
ld_1 <- paste0("chr", ld_1)
ld_1 <- gsub(":", "_", ld_1)
genotypes_ld1 <- genotypes %>% filter(grepl(paste(ld_1, collapse = "|"), het_variants))
# For the variants which don't have correlation over 0.01
# What is the maximum R2 value (for write up)
ld_distinct_pairs %>% filter(R2 < 0.01) %>% pull(R2) %>% max()
# Proportion of the variants with a higher LD than 0.01 which are carried together in the same participants
# Should all be 2 (i.e if the LD is 1 then the number of het_variants should be two as the pairs should always be carried together)
table(genotypes_ld1$n_het_variants)
pairs_over0.01 <- ld %>%
  filter(SNP_A != SNP_B & R2 > 0.01) %>% # Filter as you have done
  mutate(SNP_A_new = pmin(SNP_A, SNP_B),          # Ensure consistent ordering
         SNP_B_new = pmax(SNP_A, SNP_B)) %>% 
  select(SNP_A = SNP_A_new, SNP_B = SNP_B_new, R2) %>% # Select new columns
  distinct()  
write.table(pairs_over0.01, paste0(gene, "_corrpairs_over0.01.tsv"), sep = "\t", row.names =F, quote =F)

  # LD matrix of all those showing LD over 0.01
ld_over0.01_matrix <- ld_matrix[row.names(ld_matrix) %in% ld_over0.01$SNP_A, colnames(ld_matrix) %in% ld_over0.01$SNP_A]

png(paste0(gene, "_priority_heatmap_over0.01.png"), width = 8, height = 6, units = "in", res = 300)
pheatmap(as.matrix(ld_over0.01_matrix),
 cluster_cols = FALSE,
cluster_rows = FALSE,
 display_numbers = FALSE,
 breaks = seq(0,1, length.out = 101),
 main = paste0(gene, ": R2 > 0.01 (n =", nrow(pairs_over0.01), " distinct ", ifelse(nrow(pairs_over0.01)==1, 'pair)', 'pairs)')))
dev.off()

} else {
   print("No variants with > 0.01 correlation")
}
}

for(i in 1:nrow(fads)) {
    gene_name <- fads$hgnc_symbol[i]
    print(gene_name)
    plot_heatmaps(gene_name)
}

