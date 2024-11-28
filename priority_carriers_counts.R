### Reading in the priority genotypes ###
# Example with FEN1

# dx download /Output/genotypes/FEN1_priority_genotypes.tsv
# dx download /Output/gene_VCF_variants/variants/pri_variants/FEN1_priority_annot_chrposrefalt.txt
# dx download /Input/FADS_cluster_UKB_pVCF.tsv

### Read in the libraries
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(optparse)

# Read in the FADS cluster co-ordinates file
fads_genes <- read.table("FADS_cluster_UKB_pVCF.tsv", sep = "\t" , header = T)
# The genotypes of the prioritised variants 
genotypes <- fread("FEN1_priority_genotypes.tsv", sep = "\t", header = F)
colnames(genotypes) <- c("CHR", "POS", "REF", "ALT", "SAMPLE", "GT")
genotypes <- genotypes %>% mutate(chr_pos_ref_alt = paste0(CHR, "_", POS, "_", REF, "_", ALT))

# How many variants have been pulled out 
length(unique(genotypes$chr_pos_ref_alt)) 
# Read in the priority variant list
pri_variants <- readr::read_lines("FEN1_priority_annot_chrposrefalt.txt")

# Are there the same number variants pulled out than in the priority list 
length(pri_variants) == length(unique(genotypes$chr_pos_ref_alt))

# Pull out the extra variants which have been extracted (should have the same position as a priority variant)

unique(genotypes$chr_pos_ref_alt)[unique(genotypes$chr_pos_ref_alt) %in% pri_variants == FALSE]

# Filter to just the variants in the priority variant list 

genotypes <- genotypes %>% filter(chr_pos_ref_alt %in% pri_variants)

# Find out how many participants carry a prioritised variant 
table(genotypes$GT)
# There are some . and some ./. - how to deal with these? 
# - remove the individuals with a missing phenotype at any prioritised variant 
# - As cannot be certain if these individuals carry the variant or not 

genotypes <- genotypes %>% filter(GT != "./." & GT != ".")

# See how many carriers per variant

variant_carrier_summary <- genotypes %>% group_by(chr_pos_ref_alt) %>% 
summarise(althet_carriers = sum(GT == "0/1"),
althom_carriers = sum(GT == "1/1"),
refhom_carriers= sum(GT=="0/0"))

variant_carrier_summary <- variant_carrier_summary %>%
separate(chr_pos_ref_alt, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
mutate(POS = as.numeric(POS))

# How many variants carried per participant

participant_variant_summary <- genotypes %>% 
group_by(SAMPLE) %>% 
summarise(num_althet_vars = sum(GT=="0/1"),
num_althom_vars = sum(GT=="1/1"),
num_altref_vars = sum(GT=="0/0"))

ggplot(participant_variant_summary, 
aes(x = num_althet_vars, fill = as.factor(num_althet_vars))) + 
geom_bar(stat="count")+  
geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5) + 
theme_minimal() + 
labs(x = "Number of FEN1 prioritised variants carried per individual", y = "Number of participants") +
theme(legend.position="none")

# Extract the carriers 

carriers <- participant_variant_summary %>% 
filter(num_althom_vars > 0 | num_althet_vars > 0)

# Extract non carriers
# This excludes those who are reference homozygous in almost all variants, but have a missing genotype at one (n = 13)
noncarriers <- participant_variant_summary %>% 
filter(num_altref_vars == length(pri_variants))

# Pull out the samples with no alternate heterozygous variants or alternate homozygous variants
samples_with_no_althet_vars <- participant_variant_summary %>%
    filter(num_althet_vars == 0 & num_althom_vars == 0) %>%
    pull(SAMPLE)

# Step 2: Find which samples in samples_with_no_althet_vars are not in noncarriers$SAMPLE
samples_not_in_noncarriers <- samples_with_no_althet_vars[!samples_with_no_althet_vars %in% noncarriers$SAMPLE]

# Establish how many of reference homozygous genotypes these samples have 
# (These have less than the number of priority variants)
participant_variant_summary %>% 
filter(SAMPLE %in% samples_not_in_noncarriers) %>% 
select(SAMPLE, num_altref_vars)

## Get genotype information for the carriers

carriers_genotypes <- genotypes %>% 
filter(SAMPLE %in% carriers$SAMPLE) %>% 
group_by(SAMPLE) %>% 
summarise(n_het_variants = sum(GT != "0/0"),
n_althom_variants = sum(GT == "1/1"),
het_variants = paste(chr_pos_ref_alt[GT!= "0/0"], collapse = ':'),
hom_variants = paste(chr_pos_ref_alt[GT=="1/1"], collapse = ",")
)

# Create a  FEN1 carrier 'phenotype' for use in plotting

carriers_info <- rbind(carriers_genotypes %>%
 mutate(status="carrier") %>% 
 select(SAMPLE, status),
 noncarriers %>% 
 mutate(status = "non-carrier") %>% 
 select(SAMPLE, status))

# Save the carrier variant information, the carrier 'phenotype', the carrier summary and the variant summary

write.table(paste0(gene, "_participant_var_summary.tsv"), participant_variant_summary, sep = '\t', row.names = F, quote = F)
write.table(paste0(gene, "_variant_count_summary.tsv"), variant_carrier_summary, sep = '\t', row.names = F, quote = F)
write.table(paste0(gene, "_carrier_info.tsv"), carriers_info, sep = "\t", row.names = F, quote = F)
write.table(paste0(gene, "_carriers_genotypes.tsv"), carriers_genotypes, sep = "\t", row.names = F, quote =F)