library(ggalt)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr) # Have to install this in the post workbench
# dx download the variants 
# dx download /Output/gene_VCF_variants/variants/*QC_variants.tsv
# dx download /Input/FADS_cluster_UKB_pVCF.tsv

# Read in the FADS coordinates used to index the VCF files 
fads <- read.table("FADS_cluster_UKB_pVCF.tsv", sep = "\t", header = T)
# Read in the variants files
for (i in c(1:nrow(fads))) {
    gene <- fads$hgnc_symbol[i]
    print(gene)
    print(paste0("gene_co-ordinates: ", fads$start_position[i], "-", fads$end_position[i]))
    gene_table <- read.table(paste0(gene, "_QC_variants.tsv"), sep = "\t")
    colnames(gene_table) <- c("CHR", "POS", "REF", "ALT", "FILTER", "AC", "AN")
    gene_table <- gene_table %>% mutate(chr_pos_ref_alt = paste0(CHR, "_", POS, "_", REF, "_", ALT), 
    AF = AC/AN)
    # Remove variants with an AC == 0 (these variants had a low allele count before sample QC)
    # And the participants which carried them have been removed- hence allele count == 0
    variant_counts <- data.frame(GENE=gene, total = length(gene_table$chr_pos_ref_alt))
    gene_AC_zero <- gene_table %>% filter(AC==0) %>% pull(chr_pos_ref_alt)
    variant_counts <- variant_counts %>% mutate(AC_zero=length(unique(gene_AC_zero)))
    gene_table <- gene_table %>% filter(AC != 0)
    gene_n <- length(unique(gene_table$chr_pos_ref_alt))
    variant_counts <- variant_counts %>% mutate(AC_over0_postQC=gene_n)
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
  count(range) %>%
   mutate(GENE=paste0(gene, ' (n=', gene_n, ')'))
    AC_counts <- gene_table %>% mutate(
range = case_when(
  AC == 1 ~ "1",
  AC > 1 & AC < 5 ~ "< 5",
  AC >= 5 & AC < 10 ~ "< 10",
  AC >= 10 & AC < 50 ~ "< 50",
  AC >= 50 & AC < 100 ~ "< 100",
  AC >= 100 & AC < 1000 ~ "< 1000",
  AC >= 1000 & AC < 10000 ~ "< 10000",
  AC >= 10000 ~ ">= 10000"
)
  ) %>%
  count(range) %>% 
  mutate(GENE=paste0(gene, ' (n=', gene_n, ')'), range = factor(range, levels = c("1", "< 5", "< 10", "< 50", "< 100", "< 1000", "< 10000", ">= 10000")))

assign(paste0(gene, "_AF_counts"), frequency_counts)
assign(paste0(gene, "_AC_counts"), AC_counts)
assign(paste0(gene, "_AC_zero_postqc"), gene_AC_zero)
assign(paste0(gene, "_variant_counts"), variant_counts)
}
all_AF_counts <- rbind(MYRF_AF_counts, TMEM258_AF_counts, FEN1_AF_counts, FADS2_AF_counts,FADS1_AF_counts,  FADS3_AF_counts)
all_AC_counts <- rbind(MYRF_AC_counts, TMEM258_AC_counts, FEN1_AC_counts, FADS2_AC_counts,FADS1_AC_counts,  FADS3_AC_counts)


gene_af_plot <- ggplot(all_AF_counts, aes(x = range, y = n)) + 
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 2) + 
  labs(x = "UKB Allele Frequency Range", y = "Count", title = "UKB Allele Frequency (post-QC)") + 
  theme_minimal() +
  facet_wrap(~GENE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

gene_ac_plot <- ggplot(all_AC_counts, aes(x=range, y = n)) + 
  geom_bar(stat="identity", fill = "skyblue")+
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 2) + 
  labs(x = "UKB Allele Count", y = "Count", title="UKB Allele Count (post-QC)") + 
  theme_minimal() +
  facet_wrap(~GENE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "AF_barplots.png", gene_af_plot, width = 10, height = 8, device = "png", dpi = 300)
ggsave(filename = "AC_barplots.png",gene_ac_plot, width = 10, height = 8, device = "png", dpi = 300)

# Write out the list of variants with Allele Count == 0 in each gene 
readr::write_lines(FADS1_AC_zero_postqc, "FADS1_postQC_ACzero.txt")
readr::write_lines(FADS2_AC_zero_postqc, "FADS2_postQC_ACzero.txt")
readr::write_lines(FADS3_AC_zero_postqc, "FADS3_postQC_ACzero.txt")
readr::write_lines(FEN1_AC_zero_postqc, "FEN1_postQC_ACzero.txt")
readr::write_lines(MYRF_AC_zero_postqc, "MYRF_postQC_ACzero.txt")
readr::write_lines(TMEM258_AC_zero_postqc, "TMEM258_postQC_ACzero.txt")

AC_counts <- rbind(FADS1_variant_counts, FADS2_variant_counts, FADS3_variant_counts, FEN1_variant_counts, MYRF_variant_counts, TMEM258_variant_counts)
AC_counts_long <-  AC_counts  %>%
    pivot_longer(cols = -GENE, 
                 names_to = "Category", 
                 values_to = "Count") %>%
                 mutate(Category = factor(Category, levels = c("total", "AC_zero", "AC_over0_postQC")))

ac_zero_counts_plt <- ggplot(AC_counts_long, aes(x=GENE, y = Count, fill = Category)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.8))+
  geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.7, size = 2) +
  labs(x = "GENE", y = "Count") + 
  theme_minimal() +
  theme(legend.position = "top")+ 
  scale_fill_manual(
    values = c("total" = "gray",
    "AC_zero" = "red",
    "AC_over0_postQC" = "#C77CFF"
    ),
labels = c("total"="Total variants (post-QC)",
"AC_zero"="QC variants AC = 0",
"AC_over0_postQC"="QC variants AC > 0" ))

ggsave(filename = "AC_zero_postqc_barplot.png",ac_zero_counts_plt, width = 10, height = 8, device = "png", dpi = 300)

# dx upload *_postQC_ACzero.txt
# dx upload AC_zero_postqc_barplot.png