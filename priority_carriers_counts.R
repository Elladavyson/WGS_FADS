### Reading in the priority genotypes ###
# Example with FEN1

# dx download /Output/genotypes/FEN1_priority_genotypes.tsv
# dx download /Output/genotypes/FEN1_extra_priority_genotypes.tsv
# dx download /Output/gene_VCF_variants/variants/pri_variants/chrpos/chrposrefalt_gene_canonical/FEN1_priority_annot_chrposrefalt.txt
# dx download /Input/FADS_cluster_UKB_pVCF.tsv

### Read in the libraries
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c('-g','--gene'), type = 'character', help = 'Gene looking at', metavar = "GENE")
parser <- add_option(parser, c('-e','--extra'), type = 'character', help = 'Are there extra genotypes pulled out', metavar = "EXTRA")

opt=parse_args(parser)
print(opt)

gene = opt$gene
extra = opt$extra
sink(paste0(gene, "_summarise_carriers.log"))
print("--------------------------------------------------")
print("READING IN DATA")
print("--------------------------------------------------")
# Read in the FADS cluster co-ordinates file
print("Reading in files")
fads_genes <- read.table("FADS_cluster_UKB_pVCF.tsv", sep = "\t" , header = T)
# The genotypes of the prioritised variants 
genotypes <- fread(paste0(gene, "_priority_genotypes.tsv"), sep = "\t", header = F)
colnames(genotypes) <- c("CHR", "POS", "REF", "ALT", "SAMPLE", "GT")
genotypes <- genotypes %>% mutate(chr_pos_ref_alt = paste0(CHR, "_", POS, "_", REF, "_", ALT))
# If there are extra genotypes pulled out due to the canonical issue for the gene
# Read these in and bind the genotypes together
if(extra =="Yes"){
    print("Appending the exta genotypes")
    extra_genotypes <- fread(paste0(gene, "_extra_priority_genotypes.tsv"), sep = "\t", header =F)
    # Make sure the extra genotypes ARE extra (not creating any duplications)
    colnames(extra_genotypes) <- c("CHR", "POS", "REF", "ALT", "SAMPLE", "GT")
    extra_genotypes <- extra_genotypes %>% mutate(chr_pos_ref_alt = paste0(CHR, "_", POS, "_", REF, "_", ALT))
    extra_genotypes <- extra_genotypes %>% filter(chr_pos_ref_alt %in% unique(genotypes$chr_pos_ref_alt) == FALSE)
    genotypes <- rbind(genotypes, extra_genotypes)
}

print("--------------------------------------------------")
print("CHECKING VARIANTS PULLED OUT")
print("--------------------------------------------------")
# How many variants have been pulled out 
print("How many prioritised variants have been pulled out of BCFTOOLS using CHR and POS")
length(unique(genotypes$chr_pos_ref_alt)) 
# Read in the priority variant list
pri_variants <- readr::read_lines(paste0(gene, "_priority_annot_chrposrefalt.txt"))
print("How many prioritised variants are there in the list")
length(pri_variants)
pri_variants_df <- data.frame(chr_pos_ref_alt = pri_variants) %>% 
separate(chr_pos_ref_alt, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE)
# Are there the same number variants pulled out than in the priority list 
print("Are the number of variants pulled out by BCFTOOLS and the number of those in the list the same?")
length(pri_variants) == length(unique(genotypes$chr_pos_ref_alt))

# Pull out the extra variants which have been extracted (should have the same position as a priority variant)
print("The extra variants extracted")
extra_variants_bcftools <- genotypes %>% filter(chr_pos_ref_alt %in% pri_variants == FALSE) %>%
distinct(chr_pos_ref_alt, .keep_all = TRUE)
table(extra_variants_bcftools$POS %in% pri_variants_df$POS)
# Filter to just the variants in the priority variant list 
print("Filtering to just the variants in the variant list")
genotypes <- genotypes %>% filter(chr_pos_ref_alt %in% pri_variants)
# Check that all the priority variants are present 
if(table(pri_variants %in% genotypes$chr_pos_ref_alt)==FALSE) {
    stop("Not all priority variants pulled out, investigate")
} else {
    print('All priority variants pulled out')
}

print("--------------------------------------------------")
print("SUMMARISING GENOTYPE COUNTS PER PARTICIPANT AND PER VARIANT")
print("--------------------------------------------------")
# Find out how many participants carry a prioritised variant 
print("The distribution of genotypes in the prioritised variants")
table(genotypes$GT)
# There are some . and some ./. - how to deal with these? 
# - remove the individuals with a missing phenotype at any prioritised variant 
# - As cannot be certain if these individuals carry the variant or not 
print("Remove the missing genotypes")
genotypes <- genotypes %>% filter(GT != "./." & GT != ".")

# See how many carriers per variant
print("See how many carriers each variant has and what type")
variant_carrier_summary <- genotypes %>% group_by(chr_pos_ref_alt) %>% 
summarise(althet_carriers = sum(GT == "0/1"),
althom_carriers = sum(GT == "1/1"),
refhom_carriers= sum(GT=="0/0"))

variant_carrier_summary <- variant_carrier_summary %>%
separate(chr_pos_ref_alt, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE) %>%
mutate(POS = as.numeric(POS))

print(variant_carrier_summary)

var_sum_plt <- ggplot(variant_carrier_summary, aes(x = althet_carriers)) + 
geom_histogram() + 
labs(x = "Number of alternate heterozygous carriers per variant", y = "Number of variants", title = paste0(gene, ": prioritised variant carriers")) + 
theme_minimal() 

# How many variants carried per participant
print("How many variants carried by each individual")
participant_variant_summary <- genotypes %>% 
group_by(SAMPLE) %>% 
summarise(num_althet_vars = sum(GT=="0/1"),
num_althom_vars = sum(GT=="1/1"),
num_altref_vars = sum(GT=="0/0"))
print("How many variants carried in alternate heterozygous per person")
table(participant_variant_summary$num_althet_vars)
print("How many variants carried in alternate homozygous per person")
table(participant_variant_summary$num_althom_vars)
part_var_plt <- ggplot(participant_variant_summary %>% filter(num_althet_vars !=0), 
aes(x = as.factor(num_althet_vars), fill = as.factor(num_althet_vars))) + 
geom_bar(stat="count")+  
geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5) + 
theme_minimal() + 
labs(x = paste0("Number of ", gene, " prioritised variants carried per individual"), y = "Number of participants", title = paste0(gene, ": prioritised variant carriers")) +
theme(legend.position="none")

print("--------------------------------------------------")
print("EXTRACTING THE CARRIERS AND NON-CARRIERS")
print("--------------------------------------------------")
# Extract the carriers 
print("Extracting the carriers")
carriers <- participant_variant_summary %>% 
filter(num_althom_vars > 0 | num_althet_vars > 0)
print(paste0("There are ", length(unique(carriers$SAMPLE)), " who carry at least one of prioritised variants in ", gene))
# Extract non carriers
# This excludes those who are reference homozygous in almost all variants, but have a missing genotype at one (n = 13)
noncarriers <- participant_variant_summary %>% 
filter(num_altref_vars == length(pri_variants))
print(paste0("There are ", length(unique(noncarriers$SAMPLE)), " who do NOT carry any of prioritised variants in ", gene))
# Pull out the samples with no alternate heterozygous variants or alternate homozygous variants
samples_with_no_althet_vars <- participant_variant_summary %>%
    filter(num_althet_vars == 0 & num_althom_vars == 0) %>%
    pull(SAMPLE)

# Step 2: Find which samples in samples_with_no_althet_vars are not in noncarriers$SAMPLE
samples_not_in_noncarriers <- samples_with_no_althet_vars[!samples_with_no_althet_vars %in% noncarriers$SAMPLE]

# Establish how many of reference homozygous genotypes these samples have 
# (These have less than the number of priority variants)
print("The samples which have at least one missing genotype for a priority variant and therefore we cannot ascertain carrier status or not for sure")
print("Removing these samples")
participant_variant_summary %>% 
filter(SAMPLE %in% samples_not_in_noncarriers) %>% 
select(SAMPLE, num_altref_vars)

## Get genotype information for the carriers
print("Extracting the genotype information for carriers to save")
carriers_genotypes <- genotypes %>% 
filter(SAMPLE %in% carriers$SAMPLE) %>% 
group_by(SAMPLE) %>% 
summarise(n_het_variants = sum(GT != "0/0"),
n_althom_variants = sum(GT == "1/1"),
het_variants = paste(chr_pos_ref_alt[GT!= "0/0"], collapse = ':'),
hom_variants = paste(chr_pos_ref_alt[GT=="1/1"], collapse = ",")
)

# Create a  FEN1 carrier 'phenotype' for use in plotting
print(paste0("Creating a variable for the ", gene, " which denotes participants which carry a variant as 1 and those who do not as 0"))
carriers_info <- rbind(carriers_genotypes %>%
 mutate(status="carrier") %>% 
 select(SAMPLE, status),
 noncarriers %>% 
 mutate(status = "non-carrier") %>% 
 select(SAMPLE, status))

table(carriers_info$status)
# Establish per variant information 

# Establish the compound heterozygous variants 
poten_comp_het_variants <- carriers_genotypes %>% 
filter(n_het_variants > 1 & n_althom_variants == 0 ) %>% 
pull(het_variants) %>% strsplit(., ':') %>%
unlist() %>% 
unique()

variant_carrier_summary <- variant_carrier_summary %>% 
mutate(zygosity = case_when(chr_pos_ref_alt %in% poten_comp_het_variants ~ "Inconc_comp_het",
althet_carriers > 0 & althom_carriers > 0 ~ "Alt_hom_and_alt_het",
althet_carriers > 0 & althom_carriers ==0 ~ "Alt-het",
althet_carriers == 0 & althom_carriers > 0 ~ "Alt-hom",
TRUE ~ "No valid zygosity"))

variant <- variant_carrier_summary %>% mutate(colour_for_lolliplot = 
case_when(zygosity == "Alt-het" ~ "#90EE90",
zygosity == "Alt_hom_and_alt_het" ~ "#006400", 
zygosity == "Alt-hom" ~ "#ADD8E6",
zygosity == "Inconc_comp_het" ~  "lightpink",
TRUE ~ "No valid zygosity"))

variant_info <- variant %>% 
select(c(chr_pos_ref_alt, CHR, POS, REF, ALT, zygosity, colour_for_lolliplot, althet_carriers))

# Save the carrier variant information, the carrier 'phenotype', the carrier summary and the variant summary

print("--------------------------------------------------")
print("SAVING")
print("--------------------------------------------------")

print(paste0("Saving summary of the number of variants per participant to: ", gene, "_participant_var_summary.tsv"))
print(paste0("Saving summary of the number of carriers per variant to: ", gene, "_variant_count_summary.tsv"))
print(paste0("Saving the carrier status file to: ", gene, "_carriers_info.tsv"))
print(paste0("Saving the genotype information about the carriers to: ", gene, "_carriers_genotypes.tsv"))
print(paste0("Saving per variant information about zygosity to: ", gene, "_variant_zygosity.tsv"))

write.table(participant_variant_summary, paste0(gene, "_participant_var_summary.tsv"), sep = '\t', row.names = F, quote = F)
write.table(variant_carrier_summary, paste0(gene, "_variant_count_summary.tsv"), sep = '\t', row.names = F, quote = F)
write.table(carriers_info, paste0(gene, "_carrier_info.tsv"),  sep = "\t", row.names = F, quote = F)
write.table(carriers_genotypes, paste0(gene, "_carriers_genotypes.tsv"),  sep = "\t", row.names = F, quote =F)
write.table(variant_info, paste0(gene, "_variant_zygosity.tsv"), sep = "\t", row.names = F, quote = F)

# Save the graphs 
print(paste0("Graph summarising the number of variants per participants saved to : ", gene, "_pri_carriers_summary.png"))
print(paste0("Graph summarising the number of participants per variant saved to: ", gene, "_pri_variant_summary.png"))

ggsave(filename = paste0(gene, "_pri_carriers_summary.png"), part_var_plt, width = 10, height = 6, device = "png", dpi = 300)
ggsave(filename = paste0(gene, "_pri_variant_summary.png"), var_sum_plt, width = 10, height = 6, device = "png", dpi = 300)
