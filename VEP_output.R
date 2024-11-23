
######## LIBRARIES ##############
.libPaths(c(.libPaths(), "/home/s2112198/edavyson/x86_64-pc-linux-gnu-library/4.4"))
library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(UpSetR)
library(tidyr)
library(optparse)

######## OPTION SET UP  ##############
parse <- OptionParser()

option_list <- list(
    make_option('--gene', type = "character", help = "Gene name", action = "store")
)
args=commandArgs(trailingOnly = TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args = args)

gene <- opt$gene

file_dir <-paste0(getwd(), "/", gene, "_vepoutput_proc/")
# Check if the directory exists
if (!dir.exists(file_dir)) {
  # If it doesn't exist, create the directory
  dir.create(file_dir)
  message("Directory created for output files: ", file_dir)
} else {
  message("File Directory already exists: ", file_dir)
}

sink(paste0(file_dir, gene, "_vepoutput.log"))

  vep_consequence <- factor(c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", 
  "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "feature_elongation", 
  "feature_truncation", "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant", 
  "splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
   "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
    "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant", 
    "NMD_transcript_variant", "non_coding_transcript_variant", "coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", 
    "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification", "regulatory_region_variant", "intergenic_variant"),
      levels = c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", 
             "stop_gained", "frameshift_variant", "stop_lost", "start_lost", 
             "transcript_amplification", "feature_elongation", "feature_truncation", 
             "inframe_insertion", "inframe_deletion", "missense_variant", 
             "protein_altering_variant", "splice_donor_5th_base_variant", 
             "splice_region_variant", "splice_donor_region_variant", 
             "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant", 
             "start_retained_variant", "stop_retained_variant", "synonymous_variant", 
             "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", 
             "3_prime_UTR_variant", "non_coding_transcript_exon_variant", 
             "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant", 
             "coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", 
             "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", 
             "regulatory_region_ablation", "regulatory_region_amplification", 
             "regulatory_region_variant", "intergenic_variant"),
  ordered = TRUE)

# Read in the VEP output file 
vep_output <- fread(paste0(gene, "_QC_nosamples_vepannot.vcf_annotated_hdr.tsv"), sep = "\t")
# Check how many variants in the file 
length(unique(vep_output$chr_pos_ref_alt))
# Filter to the canonical transcript
vep_output_canonical <- vep_output %>% filter(CANONICAL == "YES")
length(unique(vep_output_canonical$chr_pos_ref_alt))
# For variants with multiple consequences, filter to the most severe 
vep_output_canonical <- vep_output_canonical %>% mutate(Consequence_factor = factor(Consequence, levels = levels(vep_consequence), ordered = TRUE))
vep_output_canonical$chr_pos_ref_alt <- gsub(':', '_', vep_output_canonical$chr_pos_ref_alt)
annotated <- vep_output_canonical %>% 
group_by(chr_pos_ref_alt) %>% 
slice_min(Consequence_factor,with_ties = FALSE) %>%
ungroup()
length(unique(annotated$chr_pos_ref_alt))

# Read in variant allele frequency 
internal_AF <- read.table(paste0("~/edavyson/WGS_FADS/QC/", gene, "_QC_variants.tsv"), sep = "\t")
 colnames(internal_AF) <- c("CHR", "POS", "REF", "ALT", "FILTER", "AC", "AN")
internal_AF <- internal_AF %>% mutate(chr_pos_ref_alt = paste0(CHR, "_", POS, "_", REF, "_", ALT), 
    AF = AC/AN)

annotated_internalAF <- merge(annotated, internal_AF, by = "chr_pos_ref_alt")
annotated_internalAF <- annotated_internalAF %>%
  mutate(
    meets_filter = ifelse(AF < 0.01 & (gnomADg_AF_nfe < 0.01 | gnomADg_AF_nfe == "-"), 
                          "Yes", "No")
  )

# Filter to rare variants 
annotated_rare <- annotated_internalAF %>% filter(AF < 0.01 & (gnomADg_AF_nfe < 0.01 | gnomADg_AF_nfe == "-"))
annotated_rare$CADD_PHRED <- as.numeric(annotated_rare$CADD_PHRED)
annotated_rare$REVEL <- as.numeric(annotated_rare$REVEL)

## Setting up plots directory 

plot_dir <-paste0(getwd(), "/", gene, "_priority_plots/")
# Check if the directory exists
if (!dir.exists(plot_dir)) {
  # If it doesn't exist, create the directory
  dir.create(plot_dir)
  message("Directory created for plots: ", plot_dir)
} else {
  message("Plot Directory already exists: ", plot_dir)
}
# How many variants meet the rare criteria
af_filter_plt <- ggplot(annotated_internalAF, aes(x = meets_filter)) + geom_bar(stat="count") + 
labs(x = "Variants with AF < 0.01 and gnomADg_AF_nfe < 0.01 or -", y = "Number of variants", title = gene) + 
theme_minimal()

  png(paste0(plot_dir, "AF_filter_", gene, ".png"), 
    width = 3000, height = 2000, res = 300, type = "cairo")
af_filter_plt
dev.off()

# What is the distribution of the allele count in UKB of the rare variants?

ac_hist <- ggplot(annotated_rare, aes(x = AC)) + geom_histogram(fill = "skyblue") + 
labs(x = "Allele Count of variant", y = "Count", title = paste0(gene, ": distribution of allele count of rare variants"))

ac_hist_zoom <- ggplot(annotated_rare, aes(x = AC)) + geom_histogram(binwidth=1, fill = "skyblue", boundary = 0) + 
labs(x = "Allele Count of variant", y="Count", title = paste0(gene, ": distribution of allele count of rare variants (0-25)")) +
xlim(0,25) + stat_bin(binwidth = 1, geom = "text", aes(label = ..count..), vjust = -0.5, boundary = 0)


  png(paste0(plot_dir, "AC_hist_rare", gene, ".png"), 
    width = 3000, height = 3000, res = 300, type = "cairo")
ggarrange(ac_hist, ac_hist_zoom, nrow = 2, ncol =1 )
dev.off()

# What are the consequences annotated to the rare variants 
csq_plot <- ggplot(annotated_rare, aes(x = reorder(Consequence, table(Consequence)[Consequence]), 
                                      fill = reorder(IMPACT, table(Consequence)[Consequence]))) +
  geom_bar(stat = 'count') + theme_bw() + 
 labs(x = "Consequence", y = "N", title = paste0(gene, " rare variant VEP Consequences"), fill = "VEP IMPACT") +
  coord_flip() + geom_text(stat = 'count', aes(label = after_stat(count)), hjust = -0.1) + 
  theme(legend.position = "top")

  png(paste0(plot_dir, "CSQ_rare", gene, ".png"), 
    width = 5000, height = 2000, res = 300, type = "cairo")
csq_plot
dev.off()

# pathogenicity filters

## LoFTee

lof_plt  <- ggplot(annotated_rare %>% filter(LoF != '-'), aes(x = LoF, fill = LoF)) + geom_bar(stat="count") +
  labs(x = "LoFTEE Score", y = "Count",title = paste0("LOFTEE rating for rare variants annotated to canonical transcript of ", gene)) +
  geom_text(stat='count',aes(label=after_stat(count),vjust=-1, color = LoF)) +
  theme_minimal() + 
    scale_fill_manual(values = c("HC" = "red", "LC"="black", "OS"="lightblue", "." = "black")) +  # Color bars
  scale_color_manual(values = c("HC" = "red", "LC"="black", "OS"="lightblue", "." = "black"))+
    guides(color = "none")
png(paste0(plot_dir, "/LoFTee_", gene, ".png"), 
    width = 3000, height = 2000, res = 300, type = "cairo")
lof_plt
dev.off()

## CADD
cadd_plt <- ggplot(annotated_rare, aes(x = CADD_PHRED)) + 
  geom_histogram(binwidth = 1) +  # Specify binwidth for better control
  labs(x = "CADD PHRED", y = "Count", 
       title = paste0("CADD scores for rare variants annotated to canonical transcript of ", gene)) + 
  theme_minimal() +
  geom_vline(xintercept = 20, linetype = "dashed", color = "red") 
png(paste0(plot_dir, "/CADD_", gene ,".png"), 
    width = 3000, height = 2000, res = 300, type = "cairo")
cadd_plt
dev.off() 

## REVEL 

png(paste0(plot_dir, "/REVEL_", gene, ".png"), 
    width = 3000, height = 2000, res = 300, type = "cairo")
ggplot(annotated_rare, aes(x = as.numeric(REVEL), fill = as.factor(REVEL >= 0.5))) + geom_histogram() +
  labs(x = "REVEL", y = "Count", title = paste0("REVEL scores for rare variants annotated to canonical transcript of ", gene)) +
  theme_minimal() +
  xlim(0,1) +
  geom_vline(xintercept=0.5, linetype = "dashed", color = 'red') 
  # 6350 have '-'
dev.off()

table(annotated_rare$CADD_PHRED > 30)
table(annotated_rare$IMPACT)
table(annotated_rare$LoF)
table(annotated$REVEL > 5)

priority <- annotated_rare %>%
  mutate(priority_criteria = case_when(
    LoF == "HC" | LoF == "OS" ~ "LoF: HC/OS",
    TRUE ~ NA_character_
  )) %>%
  mutate(priority_criteria = case_when(
    CADD_PHRED >= 20 ~ ifelse(is.na(priority_criteria), "CADD >= 20", paste(priority_criteria, "CADD > 20", sep = ", ")),
    TRUE ~ priority_criteria
  )) %>%
  mutate(priority_criteria = case_when(
    IMPACT == "HIGH" | IMPACT == "MODERATE" ~ ifelse(is.na(priority_criteria), "IMPACT: HIGH/MODERATE", paste(priority_criteria, "IMPACT: HIGH/MODERATE", sep = ", ")),
    TRUE ~ priority_criteria
  )) %>%
  mutate(priority_criteria = case_when(
    REVEL>= 0.5~ ifelse(is.na(priority_criteria), "REVEL >= 0.5", paste(priority_criteria,"REVEL >= 0.5", sep = ", ")),
    TRUE ~ priority_criteria
  )) %>% 
  filter((LoF == "HC" | LoF == "OS") |
         (CADD_PHRED >= 20) |
         (IMPACT == "HIGH" | IMPACT == "MODERATE") |
         (REVEL >= 0.5))

# Create the list of annotations for the UpSetR plot (only those with variants)
CADD_lst <- annotated_rare %>%
  filter(CADD_PHRED > 30) %>% 
  pull(chr_pos_ref_alt) %>% 
  unique()

LoF_lst <- annotated_rare %>%
  filter((LoF == "HC" | LoF == "OS")) %>% 
  pull(chr_pos_ref_alt) %>% 
  unique()

VEP_IMPACT_lst <- annotated_rare %>%
  filter((IMPACT == "HIGH" | IMPACT == "MODERATE")) %>% 
  pull(chr_pos_ref_alt) %>% 
  unique()
REVEL_lst <- annotated_rare %>%
  filter(REVEL >= 0.5) %>% 
  pull(chr_pos_ref_alt) %>% 
  unique()

VEP_annotations <- list()

if (length(CADD_lst) > 0) {
  VEP_annotations$CADD <- CADD_lst
}
if (length(LoF_lst) > 0) {
  VEP_annotations$LOFTEE <- LoF_lst
}
if (length(VEP_IMPACT_lst) > 0) {
  VEP_annotations$VEP_IMPACT <- VEP_IMPACT_lst
}
if (length(REVEL_lst) > 0) {
  VEP_annotations$REVEL <- REVEL_lst
}

upset_plt <- upset(fromList(VEP_annotations), nsets = length(VEP_annotations), set_size.show = TRUE, order.by = 'freq')
# Plot the Upset Plot
png(paste0(plot_dir, "/UpSet_", gene, ".png"), 
    width = 3200, height = 2000, res = 300, type = "cairo")
upset_plt
dev.off()


## Save the full list of priority variants with annotations 
write.table(priority, paste0(file_dir, gene, "_priority_annot.tsv"), row.names = F, quote =F, sep = '\t')

## Save a list of the priority variants (for extracting the genotypes)
priority_write <- priority %>%
  separate(chr_pos_ref_alt, into = c("CHROM", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>% 
  select(CHROM, POS)

# List of the identifiers (CHR-POS-REF-ALT)
readr::write_lines(priority$chr_pos_ref_alt, paste0(file_dir, gene,"_priority_annot_chrposrefalt.txt"))

# BCFTOOLS -R file is two columns CHR and POS

write.table(priority_write, paste0(file_dir, gene, "_priority_annot_chrpos.tsv"),
                                  col.names = F, row.names = F, quote =F, sep = '\t')

sink()