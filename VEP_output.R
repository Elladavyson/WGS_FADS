
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
library(grid)
library(gridExtra)

######## OPTION SET UP  ##############
parse <- OptionParser()

option_list <- list(
    make_option('--gene', type = "character", help = "Gene name", action = "store")
)
args=commandArgs(trailingOnly = TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args = args)

gene <- opt$gene

file_dir <-paste0(getwd(), "/AC_overzero_06_01/", gene, "_vepoutput_proc/")
# Check if the directory exists
if (!dir.exists(file_dir)) {
  # If it doesn't exist, create the directory
  dir.create(file_dir)
  message("Directory created for output files: ", file_dir)
} else {
  message("File Directory already exists: ", file_dir)
}

############# LOGGING #############
sink(paste0(file_dir, gene, "_vepoutput_", Sys.Date(), ".log"))
print(paste0("Results from this script saved to: ", file_dir))
############ READING IN DATA #############
####### ALL VARIANTS (POST-SAMPLE AND VARIANT QC)
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
print("Reading in the VEP output file")
# Read in the VEP output file 
vep_output <- fread(paste0("VEP_output/", gene, "_QC_nosamples_vepannot.vcf_annotated_hdr.tsv"), sep = "\t")
print("How many variant annotations back from VEP")
# Check how many variants in the file 
length(unique(vep_output$chr_pos_ref_alt))
# Read in variant allele frequency 
print(paste0("Reading in the internal Allele Frequency file from ", "~/edavyson/WGS_FADS/QC/", gene, "_QC_variants.tsv"))
internal_AF <- read.table(paste0("~/edavyson/WGS_FADS/QC/", gene, "_QC_variants.tsv"), sep = "\t")
 colnames(internal_AF) <- c("CHR", "POS", "REF", "ALT", "FILTER", "AC", "AN")
internal_AF <- internal_AF %>% mutate(chr_pos_ref_alt = paste0(CHR, "_", POS, "_", REF, "_", ALT), 
    AF = AC/AN)
######### REMOVING VARIANTS WITH AC == 0 ###########
print(paste0("Removing the variants which have an allele count of 0 after variant and sample QC, n = ", internal_AF %>% 
filter(AC==0) %>% 
pull(chr_pos_ref_alt) %>% 
unique() %>% length()))

vep_output <- vep_output %>% 
mutate(chr_pos_ref_alt = gsub(":", "_", chr_pos_ref_alt)) %>%
filter(chr_pos_ref_alt %in% (internal_AF %>%
filter(AC !=0) %>% pull(chr_pos_ref_alt)))
print(paste0("Returned variants from VEP with no AC ==0: ", length(unique(vep_output$chr_pos_ref_alt))))

print("Filtering to the canonical transcript")
# Filter to the canonical transcript
vep_output_canonical <- vep_output %>% filter(CANONICAL == "YES" & SYMBOL == gene)
length(unique(vep_output_canonical$chr_pos_ref_alt))
# For variants with multiple consequences, filter to the most severe 
print("If multiple consequences per variant, filtering to the most severe")
vep_output_canonical <- vep_output_canonical %>% mutate(Consequence_factor = factor(Consequence, levels = levels(vep_consequence), ordered = TRUE))
vep_output_canonical$chr_pos_ref_alt <- gsub(':', '_', vep_output_canonical$chr_pos_ref_alt)

annotated <- vep_output_canonical %>% 
group_by(chr_pos_ref_alt) %>% 
slice_min(Consequence_factor,with_ties = FALSE) %>%
ungroup()
length(unique(annotated$chr_pos_ref_alt))

print("Merging the VEP output with the UKB allele frequency file")
annotated_internalAF <- merge(annotated, internal_AF, by = "chr_pos_ref_alt")
annotated_internalAF <- annotated_internalAF %>%
  mutate(
    meets_filter = ifelse(AF < 0.01 & (gnomADg_AF_nfe < 0.01 | gnomADg_AF_nfe == "-"), 
                          "Yes", "No")
  )

# Filter to rare variants 
print("Filtering to rare variants based off UKB AF < 0.01 AND gnomAD AF < 0.01 (or absent in gnomAD)")
annotated_rare <- annotated_internalAF %>% filter(AF < 0.01 & (gnomADg_AF_nfe < 0.01 | gnomADg_AF_nfe == "-"))
# Converting the CADD and REVEL scores to numerical variables
annotated_rare$CADD_PHRED <- as.numeric(annotated_rare$CADD_PHRED)
annotated_rare$REVEL <- as.numeric(annotated_rare$REVEL)

########## PLOTTING ############

## Setting up plots directory 
plot_dir <-paste0(getwd(), "/AC_overzero_06_01/", gene, "_priority_plots/")
# Check if the directory exists
if (!dir.exists(plot_dir)) {
  # If it doesn't exist, create the directory
  dir.create(plot_dir)
  message("Directory created for plots: ", plot_dir)
} else {
  message("Plot Directory already exists: ", plot_dir)
}
# How many variants meet the rare criteria
print(paste0("How many rare variants: ", length(unique(annotated_internalAF %>% filter(meets_filter=="Yes") %>% pull(chr_pos_ref_alt)))))
af_filter_plt <- ggplot(annotated_internalAF, aes(x = meets_filter)) + geom_bar(stat="count") + 
labs(x = "Variants with AF < 0.01 and gnomADg_AF_nfe < 0.01 or -", y = "Number of variants", title = gene) + 
theme_minimal()

png(paste0(plot_dir, "AF_filter_", gene, ".png"), 
    width = 3000, height = 2000, res = 300, type = "cairo")
af_filter_plt
dev.off()

# What is the distribution of the allele count in UKB of the rare variants?
print('The distribution of Allele Count in the rare variants (AC > 0)')
table(annotated_rare$AC)
AC_counts <- annotated_rare %>% mutate(
range = case_when(
  AC == 1 ~ "1",
  AC > 1 & AC < 5 ~ "< 5",
  AC >= 5 & AC < 10 ~ "< 10",
  AC >= 10 & AC < 50 ~ "< 50",
  AC >= 50 & AC < 100 ~ "< 100",
  AC >= 100 & AC < 1000 ~ "< 1000",
  AC >= 1000 & AC < 10000 ~ "> 1000"
)
  ) %>%
  count(range) %>% 
  mutate(range = factor(range, levels = c("1", "< 5", "< 10", "< 50", "< 100", "< 1000", "> 1000")))

AC_counts_barplot <- ggplot(AC_counts, aes(x=range, y = n)) + 
  geom_bar(stat="identity", fill = "skyblue")+
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 4) + 
  labs(x = "UKB Allele Count", y = "Count", title="UKB Allele Count (post-QC + rare variants)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

  png(paste0(plot_dir, "AC_hist_rare", gene, ".png"), 
    width = 3000, height = 3000, res = 300, type = "cairo")
AC_counts_barplot
dev.off()

# What are the consequences annotated to the rare variants 
custom_colors <- c("HIGH" = "red", "MODERATE" = "orange", "LOW" = "#90EE90", "MODIFIER" ="lightpink")
csq_plot <- ggplot(annotated_rare, aes(x = reorder(Consequence, table(Consequence)[Consequence]), 
                                      fill = reorder(IMPACT, table(Consequence)[Consequence]))) +
  geom_bar(stat = 'count') + theme_bw() + 
 labs(x = "Consequence", y = "N", title = paste0(gene, " rare variant VEP Consequences"), fill = "VEP IMPACT") +
  coord_flip() + geom_text(stat = 'count', aes(label = after_stat(count)), hjust = -0.1) + 
  theme(legend.position = "top") + scale_fill_manual(values = custom_colors)


print(paste0("Saving plot of VEP consequences to: ", plot_dir, "CSQ_rare", gene, ".png"))
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
print(paste0("Saving plot of LoFTee annotations to: ", plot_dir, "LoFTee_", gene, ".png"))
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
print(paste0("Saving plot of the CADD score distribution to: ", plot_dir, "CADD_", gene ,".png"))
png(paste0(plot_dir, "/CADD_", gene ,".png"), 
    width = 3000, height = 2000, res = 300, type = "cairo")
cadd_plt
dev.off() 

## REVEL 
print(paste0("Saving plot of the REVEL score distribution to: ", plot_dir, "REVEL_", gene, ".png"))
png(paste0(plot_dir, "/REVEL_", gene, ".png"), 
    width = 3000, height = 2000, res = 300, type = "cairo")
ggplot(annotated_rare, aes(x = as.numeric(REVEL), fill = as.factor(REVEL >= 0.5))) + geom_histogram() +
  labs(x = "REVEL", y = "Count", title = paste0("REVEL scores for rare variants annotated to canonical transcript of ", gene)) +
  theme_minimal() +
  xlim(0,1) +
  geom_vline(xintercept=0.5, linetype = "dashed", color = 'red') 
  # 6350 have '-'
dev.off()

print("How many variants have a scaled CADD score > 20")
table(annotated_rare$CADD_PHRED > 20)
print("Distribution of VEP IMPACT scores")
table(annotated_rare$IMPACT)
print("Distribution of LOFTEE High and Low Confidence")
table(annotated_rare$LoF)
print("How many variants have a REVEL score > 0.5")
table(annotated$REVEL > 0.5)

############## PRIORITISED VARIANTS ##############

print("Establishing a list of priority variants, i.e those which fit any one of the above criteria")
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

########### READ IN THE LIST OF VARAINTS WHICH ARE PRESENT IN EACH COHORT (MDD) and METABOLITE

mdd_variants <- read.table(paste0("mdd_metabol_cohortvars/", gene, "_priority_var_mdd.tsv"), sep = "\t", header =T)
metabol_variants <- read.table(paste0("mdd_metabol_cohortvars/", gene, "_priority_var_metabol.tsv"), sep = '\t', header = T)
priority_mdd <- priority %>% filter(chr_pos_ref_alt %in% mdd_variants$chr_pos_ref_alt)
priority_metabol <- priority %>% filter(chr_pos_ref_alt %in% metabol_variants$chr_pos_ref_alt)

print(paste0("There are ", length(unique(priority_mdd$chr_pos_ref_alt)), " variants prioritised with at least one criteria for loss of function present in the MDD cohort"))
print(paste0("There are ", length(unique(priority_metabol$chr_pos_ref_alt)), " variants prioritised with at least one criteria for loss of function present in the Metabolite cohort"))

# What are the consequences annotated to the prioritised variants 
csq_pri_plot <- ggplot(priority, aes(x = reorder(Consequence, table(Consequence)[Consequence]), 
                                      fill = reorder(IMPACT, table(Consequence)[Consequence]))) +
  geom_bar(stat = 'count') + theme_bw() + 
 labs(x = "Consequence", y = "N", title = paste0(gene, " prioritised variant VEP consequences"), fill = "VEP IMPACT") +
  coord_flip() + geom_text(stat = 'count', aes(label = after_stat(count)), hjust = -0.1) + 
  theme(legend.position = "top") + scale_fill_manual(values = custom_colors)
print(paste0("Saving plot of VEP consequences to: ", plot_dir, "CSQ_priority", gene, ".png"))
  png(paste0(plot_dir, "CSQ_priority", gene, ".png"), 
    width = 5000, height = 2000, res = 300, type = "cairo")
csq_pri_plot
dev.off()

csq_summary <- function(priority_dataset, cohortlab) {
  priority_dataset <- priority_dataset %>%
  group_by(Consequence) %>%
  summarise(csq_count = n(),
  IMPACT = unique(IMPACT)) %>% 
  mutate(cohort=cohortlab,
  GENE=gene)
  write.table(priority_dataset, paste0(file_dir, gene, "_priority_consequence_", cohortlab, "_cohort.tsv"), sep = "\t",
  row.names = F, quote = F)
return(priority_dataset)
}

csq_summary(priority_mdd, "MDD")
csq_summary(priority_metabol, "Metabolite")

# Create the list of annotations for the UpSetR plot (only those with variants)

upset_plot <- function(priority_dataset, cohortlab) {
CADD_lst <- priority_dataset %>%
  filter(CADD_PHRED > 30) %>% 
  pull(chr_pos_ref_alt) %>% 
  unique()

LoF_lst <- priority_dataset %>%
  filter((LoF == "HC" | LoF == "OS")) %>% 
  pull(chr_pos_ref_alt) %>% 
  unique()

VEP_IMPACT_lst <- priority_dataset %>%
  filter((IMPACT == "HIGH" | IMPACT == "MODERATE")) %>% 
  pull(chr_pos_ref_alt) %>% 
  unique()
REVEL_lst <- priority_dataset %>%
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
print(paste0("Creating an Upset plot of the priority criteria, saved: ", plot_dir, "/UpSet_", gene, "_", cohortlab, ".png"))
# Plot the Upset Plot
png(paste0(plot_dir, "/UpSet_", gene, "_", cohortlab, ".png"), width = 3200, height = 2000, res = 300, type = "cairo")
print(upset(fromList(VEP_annotations), nsets = length(VEP_annotations), set_size.show = TRUE, order.by = 'freq'))
grid.text(paste0(gene, ": Prioritised variants in the ", cohortlab, "-cohort "), x = 0.6, y = 0.97, gp = gpar(fontsize = 15))
dev.off()
}

upset_plot(priority_mdd, "MDD")
upset_plot(priority_metabol, "Metabolite")

# What is the distribution of the allele count in UKB of the prioritised variants
print('The distribution of Allele Count in the prioritised variants in the MDD cohort')
table(priority_mdd$AC)
table(priority_metabol$AC)

AC_counts <- function(priority_dataset, cohortlab) {
priority_AC_counts <- priority_dataset %>% mutate(
range = case_when(
  AC == 1 ~ "1",
  AC > 1 & AC < 5 ~ "< 5",
  AC >= 5 & AC < 10 ~ "< 10",
  AC >= 10 & AC < 50 ~ "< 50",
  AC >= 50 & AC < 100 ~ "< 100",
  AC >= 100 & AC < 1000 ~ "< 1000",
  AC >= 1000 & AC < 10000 ~ "> 1000"
)
  ) %>%
  count(range) %>% 
  mutate(range = factor(range, levels = c("1", "< 5", "< 10", "< 50", "< 100", "< 1000", "> 1000")),
  cohort = cohortlab,
  GENE = gene)
write.table(priority_AC_counts, paste0(file_dir, gene, "_priorityACcounts_", cohortlab, "_cohort.tsv"), sep = "\t",
row.names=F, quote = F)
return(priority_AC_counts)
}

AC_counts(priority_mdd, "MDD")
AC_counts(priority_metabol, "Metabolite")

priority_AC_counts <- priority %>% mutate(
range = case_when(
  AC == 1 ~ "1",
  AC > 1 & AC < 5 ~ "< 5",
  AC >= 5 & AC < 10 ~ "< 10",
  AC >= 10 & AC < 50 ~ "< 50",
  AC >= 50 & AC < 100 ~ "< 100",
  AC >= 100 & AC < 1000 ~ "< 1000",
  AC >= 1000 & AC < 10000 ~ "> 1000"
)
  ) %>%
  count(range) %>% 
  mutate(range = factor(range, levels = c("1", "< 5", "< 10", "< 50", "< 100", "< 1000", "> 1000")))

priority_AC_counts_barplot <- ggplot(priority_AC_counts, aes(x=range, y = n)) + 
  geom_bar(stat="identity", fill = "skyblue")+
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 4) + 
  labs(x = "UKB Allele Count", y = "Count", title="UKB Allele Count (prioritised variants)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

  png(paste0(plot_dir, "AC_hist_prioritised", gene, ".png"), 
    width = 3000, height = 3000, res = 300, type = "cairo")
priority_AC_counts_barplot
dev.off()

############## SAVING ###############

## Save the full list of priority variants with annotations 
print(paste0("Saving the priority variants with annotation information to: ", file_dir, gene, "_priority_annot.tsv"))
write.table(priority, paste0(file_dir, gene, "_priority_annot.tsv"), row.names = F, quote =F, sep = '\t')

## Save a list of the priority variants (for extracting the genotypes)
print(paste0("Saving a list of the priority variants CHR and POS for extracting the genotypes to: ", file_dir, gene, "_priority_annot_chrpos.tsv"))
priority_write <- priority %>%
  separate(chr_pos_ref_alt, into = c("CHROM", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>% 
  select(CHROM, POS)

# BCFTOOLS -R file is two columns CHR and POS

write.table(priority_write, paste0(file_dir, gene, "_priority_annot_chrpos.tsv"),
                                  col.names = F, row.names = F, quote =F, sep = '\t')

# List of the identifiers (CHR-POS-REF-ALT)
print(paste0("Saving a list of chr_pos_ref_alt identifiers of the priority variants to: ", file_dir, gene,"_priority_annot_chrposrefalt.txt"))
print("This is to make sure that the genotypes extracted are for the correct priority variant and not another variant with the same position but different REF and ALT, as bcftools uses just the CHR and POS variables")
readr::write_lines(priority$chr_pos_ref_alt, paste0(file_dir, gene,"_priority_annot_chrposrefalt.txt"))

print("All done! (for now :))")

sink()