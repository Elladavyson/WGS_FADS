## VEP output
library(dplyr)
library(readr)
library(readr)

# Setting up vectors 
vep_colnames <- c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence", 
  "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", 
  "Existing_variation", "IMPACT", "DISTANCE", "STRAND", "FLAGS", "CANONICAL", 
  "SOURCE", "LoF", "LoF_filter", "LoF_flags", "LoF_info", "gnomAD_100x_cov", 
  "gnomAD_10x_cov", "gnomAD_15x_cov", "gnomAD_1x_cov", "gnomAD_20x_cov", "gnomAD_25x_cov", 
  "gnomAD_30x_cov", "gnomAD_50x_cov", "gnomAD_5x_cov", "gnomAD_mean_cov", 
  "gnomAD_median_approx_cov", "gnomAD_total_DP_cov", "CADD_PHRED", "CADD_RAW", 
  "gnomADg", "gnomADg_AF_afr", "gnomADg_AF_amr", "gnomADg_AF_asj", "gnomADg_AF_eas", 
  "gnomADg_AF_fin", "gnomADg_AF_nfe", "gnomADg_AF_oth")

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

gene = FADS1
# Read in the VEP output file 
vep_output <- read.table(paste0(gene, "_QC_nosamples_vepannot.tsv"), sep = "\t")
colnames(vep_output) <- vep_colnames
# Check how many variants in the file 
length(unique(vep_output$Uploaded_variation))
# Filter to the canonical transcript
vep_output_canonical <- vep_output %>% filter(CANONICAL == "YES")
length(unique(vep_output_canonical$Uploaded_variation))
# For variants with multiple consequences, filter to the most severe 
vep_output_canonical <- vep_output_canonical %>% mutate(Consequence_factor = factor(Consequence, levels = levels(vep_consequence), ordered = TRUE))
vep_output_canonical$Uploaded_variation <- gsub('/', '_', vep_output_canonical$Uploaded_variation, )
annotated <- vep_output_canonical %>% 
group_by(Uploaded_variation) %>% 
slice_min(Consequence_factor,with_ties = FALSE) %>%
ungroup()

# Read in variant allele frequency 
internal_AF <- read.table("~/edavyson/WGS_FADS/QC/FADS1_QC_variants.tsv", sep = "\t")
 colnames(internal_AF) <- c("CHR", "POS", "REF", "ALT", "FILTER", "AC", "AN")
internal_AF <- internal_AF %>% mutate(chr_pos_ref_alt = paste0(CHR, "_", POS, "_", REF, "_", ALT), 
    AF = AC/AN)

# Filter to rare variants 
annotated_rare <- annotated %>% filter(gnomADg_AF_nfe < 0.01 | gnomADg_AF_nfe == "-")
# Get distribution of potential loss of function categories
plot_dir <-paste0(getwd(), "/", gene, "_priority_plots/")
print(paste0('Plots saved to ', plot_dir))
#dir.create(plot_dir)
csq_plot <- ggplot(annotated_rare, aes(x = reorder(Consequence, table(Consequence)[Consequence]), 
                                      fill = reorder(IMPACT, table(Consequence)[Consequence]))) +
  geom_bar(stat = 'count') + theme_bw() + 
 labs(x = "Consequence", y = "N", title = paste0(gene, " VEP Consequences")) +
  coord_flip() + geom_text(stat = 'count', aes(label = after_stat(count)), hjust = -0.1)

  png(paste0(plot_dir, "CSQ_rare", gene, ".png"), 
    width = 3000, height = 2000, res = 300, type = "cairo")
csq_plot
dev.off()

# pathogenicity filters

table(annotated_rare$CADD_PHRED > 30)
table(annotated_rare$IMPACT)
table(annotated_rare$LoF)