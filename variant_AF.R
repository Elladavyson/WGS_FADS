library(ggalt)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)

# dx download the variants 
# dx download /Output/gene_VCF_variants/variants/*

# Read in the FADS coordinates used to index the VCF files 
fads <- read.table("FADS_cluster_UKB_pVCF.tsv", sep = "\t", header = T)
# Filter to the genes which we have available now (not FADS1 or FADS2 due to the 3090 file)
fads <- fads %>% filter(hgnc_symbol %in% c("FADS1", "FADS2")==FALSE)
# Read in the variants files
for (i in c(1:nrow(fads))) {
    gene <- fads$hgnc_symbol[i]
    print(gene)
    print(paste0("gene_co-ordinates: ", fads$start_position[i], "-", fads$end_position[i]))
    gene_table <- read.table(paste0(gene, "_variants.tsv"), sep = "\t")
    colnames(gene_table) <- c("CHR", "POS", "REF", "ALT", "FILTER", "AC", "AN")
    gene_table <- gene_table %>% mutate(chr_pos_ref_alt = paste0(CHR, "_", POS, "_", REF, "_", ALT), 
    AF = AC/AN)
    print(paste0("First variant:", gene_table$chr_pos_ref_alt[1], "\n Last variant:", tail(gene_table$chr_pos_ref_alt,1)))
    assign(paste0(gene, "_variants"), gene_table)
    frequency_counts <- gene_table %>% mutate(
    range = case_when(
      AF < 0.000001 ~ "< 0.000001",
      AF < 0.00001 ~ "< 0.00001",
      AF < 0.0001 ~ "< 0.0001",
      AF < 0.001 ~ "< 0.001",
      AF < 0.01 ~ "< 0.01",
      AF < 0.1 ~ "< 0.1",
      AF < 0.5 ~ "< 0.5",
      TRUE ~ ">= 0.5"
    )
  ) %>%
  count(range)
  gene_af_plot <- ggplot(frequency_counts, aes(x=range, y = n)) + 
geom_bar(stat="identity", fill = "skyblue") +
geom_text(aes(label = n), vjust = -0.5, color = "black", size = 2) + 
labs(x = "UKB Allele Frequency Range", y = "Count", title=paste0(gene, ":", nrow(gene_table), " variants")) + theme_minimal()
assign(paste0(gene, "_AF_plt"), gene_af_plot)
}

ggsave(filename = "AF_barplots.png", ggarrange(MYRF_AF_plt, FEN1_AF_plt, FADS3_AF_plt, TMEM258_AF_plt), width = 10, height = 8, device = "png", dpi = 300)