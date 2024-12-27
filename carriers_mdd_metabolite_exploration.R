library(data.table)
library(dplyr)
library(bestNormalize) # Have to install
library(ukbnmr) # Have to install
library(purrr)
library(rlang)
library(tidyr)
library(ggplot2)
library(Hmisc) # Have to install
library(ggpubr) # Have to install 
library(ggsignif) # Have to install 

# dx download /Output/genotypes/genotype_summary/all_priority/*_carrier_info.tsv
# dx download /Output/genotypes/genotype_summary/all_priority/*_participant_var_summary.tsv
# dx download /Output/genotypes/genotype_summary/all_priority/*_variant_zygosity.tsv
# dx download /Output/genotypes/genotype_summary/all_priority/*_carriers_genotypes.tsv
# dx download /Input/FADS_cluster_UKB_pVCF.tsv
# dx download /Input/MajorDepression.ukb24262.2021-07.txt
# dx download data_participant_all.csv
# dx download data_participant_metabolomics.csv
# dx download /Output/regenie/ukb_unrel_eur_covars.covar
# dx download /Output/regenie/ukb_unrel_eur_metabol.pheno
# dx download /Output/regenie/ukb_unrel_eur_pgc3_mdd.pheno
# dx download /Output/regenie/input/annotations_FADS.tsv
# dx download /Output/regenie/input/masks_FADS.txt
# dx download /Output/regenie/input/aaf_FADS.tsv

# Looking at explicitly FADS1 Mask4 0.01 and FADS1 Mask 5 0.01
# Also look at TMEM258 for Mask4 Singletons and Mask5 singletons to see the differences in distributions

# carriers info 
print("--------------------------------------------------")
print("READING IN DATA")
print("--------------------------------------------------")
# Reading in the mask file
masks <- read.table("masks_FADS.txt", header = F)
# the annotation file
annotations <- read.table("annotations_FADS.tsv", sep = "\t", header = F)
colnames(annotations) <- c("variant", "gene", "annotation")
# Make variant variable to match the variant identifier in the VEP output
annotations <- annotations %>% 
mutate(variant_match= paste0("chr", variant)) %>% 
mutate(variant_match = gsub(":", "_", variant_match))
aaf_file <- read.table("aaf_FADS.tsv", sep = "\t", header = F)
colnames(aaf_file) <- c("variant","AF", "singleton")
# Read in the FADS cluster co-ordinates file
print("Reading in carrier information")
fads_genes <- read.table("FADS_cluster_UKB_pVCF.tsv", sep = "\t" , header = T)
for (i in 1:nrow(fads_genes)){
    gene <- fads_genes$hgnc_symbol[i]
    print(gene)
    gene_carrier <- read.table(paste0(gene, "_carrier_info.tsv"), sep = "\t", header = T)
    colnames(gene_carrier) <- c("SAMPLE", "status")
    gene_carrier <- gene_carrier %>% mutate(GENE = gene)
    zygosity <- fread(paste0(gene, "_variant_zygosity.tsv"), sep = "\t") %>% as.data.frame()
    participant_summary <- fread(paste0(gene, "_participant_var_summary.tsv"), sep="\t") %>% as.data.frame()
    participant_genotype <- fread(paste0(gene, "_carriers_genotypes.tsv"), sep = "\t", header = T) %>% as.data.frame()
    # Create two separate rows for compound heterozgyous carriers
    participant_genotype <- participant_genotype %>% 
    separate_rows(., het_variants, sep = ":")
    # Get annotations for the variants in the specific gene
    annotations_gene <- annotations %>% filter(gene == gene)
    # Merge variant annotation information with participant genotypes (creates duplicates due to compound heterozygous rows)
    participant_genotype <- left_join(participant_genotype, annotations_gene, by=c("het_variants"= "variant_match"))
    participant_genotype <- left_join(participant_genotype, aaf_file, by = c("variant"))
    # Creating per-participant variables for each mask (1 if carrier of a variant in the mask, and 0 otherwise)
    participant_genotype <- participant_genotype %>% 
    mutate(Mask1.0.01_status = ifelse(annotation == "LoF", "carrier", "non-carrier"),
    Mask2.0.01_status = ifelse(annotation == "LoF" | annotation == "H_IMPACT", "carrier", "non-carrier"),
    Mask3.0.01_status = ifelse(annotation == "LoF" | annotation == "H_IMPACT" | annotation == "REVEL", "carrier", "non-carrier"),
    Mask4.0.01_status = ifelse(annotation == "LoF" | annotation == "H_IMPACT" | annotation == "REVEL" | annotation == "CADD", "carrier", "non-carrier"),
    Mask5.0.01_status = ifelse(annotation == "LoF" | annotation == "H_IMPACT" | annotation == "REVEL" | annotation == "CADD" | annotation == "MOD_IMPACT", "carrier", "non-carrier")
    ) %>% 
    mutate(Mask1.singleton_status = ifelse(Mask1.0.01_status == "carrier" & singleton == 1, "carrier", "non_carrier"),
Mask2.singleton_status = ifelse(Mask2.0.01_status == "carrier" & singleton == 1, "carrier", "non_carrier"),
Mask3.singleton_status = ifelse(Mask3.0.01_status == "carrier" & singleton == 1, "carrier", "non_carrier"),
Mask4.singleton_status = ifelse(Mask4.0.01_status == "carrier" & singleton == 1, "carrier", "non_carrier"),
Mask5.singleton_status = ifelse(Mask5.0.01_status == "carrier" & singleton == 1, "carrier", "non_carrier"),
    )
    participant_geno_all <- participant_genotype %>% select(SAMPLE, starts_with("Mask")) 
   # For carriers of multiple varaints, ensure no duplications and they are marked as a carrier of the most deleterious variant (i.e the smallest mask)
   # Get the participants with duplicate rows (compound heterozygous)
   participant_geno_duplicated <- participant_geno_all  %>%
  filter(duplicated(SAMPLE) | duplicated(SAMPLE, fromLast = TRUE))  
  # Convert to long format, order by Mask (1-5) and assign an ID column for selection
  dup_long <- participant_geno_duplicated %>% 
    pivot_longer(cols=starts_with("Mask"),
    names_to = "Mask",
    values_to = "Mask_status") %>% 
    group_by(SAMPLE, Mask) %>% 
    arrange(SAMPLE, Mask) %>% 
    mutate(id = row_number()) %>% 
    ungroup()
    # Select the row which has the carrier status in the smallest mask
    dup_unique_id <- dup_long %>% 
    group_by(SAMPLE) %>% 
    filter(Mask_status == "carrier") %>% 
    slice_min(Mask, with_ties = TRUE) %>% 
    slice_sample(n=1) %>%
    ungroup() 
    # Merge the identified rows with the original duplicated dataframe by the ID column, to select just the rows with the smallest mask
uniq_long <- merge(dup_long, dup_unique_id %>% select(SAMPLE, id), by= c("SAMPLE", "id"))
uniq_wide <- uniq_long %>%   pivot_wider(names_from = "Mask", values_from = "Mask_status")
# Select those non duplicated to rbind to the unique rows
geno_non_duplicated <- participant_geno_all %>%
  filter(!duplicated(SAMPLE) & !duplicated(SAMPLE, fromLast = TRUE))
  # Merge the unique rows from the originally duplicated rows to the unduplicated rows
  participant_geno_merge <- bind_rows(uniq_wide, geno_non_duplicated)
# Filter to the carriers of any prioritised variant 
carriers <- gene_carrier %>% filter(status == "carrier")
# Merge with the genotype and Mask information for the carriers 
carriers <- carriers %>% 
    left_join(participant_geno_merge %>% select(-id), by = "SAMPLE")
# Establish non-carriers and define them for each mask
non_carriers <- gene_carrier %>% filter(status == "non-carrier")
non_carriers <- non_carriers %>% 
mutate(Mask1.0.01_status = "non-carrier",
    Mask2.0.01_status = "non-carrier",
    Mask3.0.01_status = "non-carrier",
    Mask4.0.01_status = "non-carrier",
    Mask5.0.01_status = "non-carrier",
    Mask1.singleton_status = "non-carrier",
    Mask2.singleton_status = "non-carrier",
    Mask3.singleton_status = "non-carrier",
    Mask4.singleton_status = "non-carrier",
    Mask5.singleton_status = "non-carrier",
    ) 
# Bind together non carriers and carriers 
gene_carrier <- rbind(carriers %>% select(SAMPLE, status, GENE, Mask1.0.01_status, Mask2.0.01_status, Mask3.0.01_status, Mask4.0.01_status, Mask5.0.01_status,
Mask1.singleton_status, Mask2.singleton_status, Mask3.singleton_status, Mask4.singleton_status, Mask5.singleton_status), non_carriers)
participant_summary <- merge(participant_summary, gene_carrier, by = "SAMPLE")
    assign(paste0(gene, "_carriers"), gene_carrier)
    assign(paste0(gene, "_variants"), zygosity)
    assign(paste0(gene, "_part_summary"), participant_summary)
    assign(paste0(gene, "_part_genotype"), participant_geno_merge)
}

print("Reading in MDD phenotype file")
mdd <- read.table("ukb_unrel_eur_pgc3_mdd.pheno", header = T)
print("Reading in the Metabolite phenotype file")
metabol <- read.table("ukb_unrel_eur_metabol.pheno", header = T)
print("The covariate file")
covars <- read.table("ukb_unrel_eur_covars.covar", header = T)
# Read in all cohort info 
all <- read.csv("data_participant_all.csv")
# colnames "eid" "p31"  "p34" "p22189" "p54_i0" "p21000_i0" "p21003_i0" "p21001_i0" "p20116_i0" "p6138_i0"  "p23449_i0" "p23444_i0" "p23450_i0" "p23457_i0" "p23459_i0" "p23443_i0"
colnames(all) <- c("f.eid", "Sex", "yob", "TDI", "AC", "ethnicity", "Age","BMI", "smoking_stat", "qualifications", "f.23449.0.0", "f.23444.0.0", "f.23450.0.0", "f.23457.0.0", "f.23459.0.0", "f.23443.0.0")


###################################################################################

# Number of prioritised variants and carriers per Mask and per gene

###################################################################################

part_summary <- function(df, carriers_df, gene) {
    df <- df %>% 
    mutate(carrier_type = case_when(num_althet_vars == 0 & num_althom_vars == 0 ~ "Non-carrier",
    num_althet_vars == 1 & num_althom_vars == 0 ~ "Alternate heterozgyous",
    num_althet_vars > 1 & num_althom_vars == 0 ~ "Compound heterozygous",
    num_althet_vars == 0 & num_althom_vars > 0 ~ "Alternate homozygous",
    num_althet_vars > 0 & num_althom_vars > 0 ~ "Alternate heterozygous & alternate homozygous"))
    return(df)
}

part_summary_list <- list(part_summary(FADS1_part_summary ,"FADS1"), part_summary(FADS2_part_summary,"FADS2"),
part_summary(FADS3_part_summary, "FADS3"), part_summary(FEN1_part_summary, "FEN1"), part_summary(MYRF_part_summary, "MYRF"),
part_summary(TMEM258_part_summary, "TMEM258"))
variant_zygosity_list <- list(FADS1_variants, FADS2_variants, FADS3_variants, FEN1_variants, MYRF_variants, TMEM258_variants)
gene_carriers_list <- list(FADS1_carriers, FADS2_carriers, FADS3_carriers, FEN1_carriers, MYRF_carriers, TMEM258_carriers)
gene_carriers_all <- do.call(rbind, gene_carriers_list)
part_summary_all <- do.call(rbind, part_summary_list)
variant_zy_all <- do.call(rbind, variant_zygosity_list)

# Participant summary information (what type of carriers etc)
part_summary_all_long <- part_summary_all %>% 
pivot_longer(cols = starts_with('Mask'),
names_to="Mask",
values_to = "Carrier_status") %>%
mutate(Mask = gsub("_status", "", Mask))

part_carriers_summary <- part_summary_all_long %>% group_by(
  GENE, Mask, Carrier_status, carrier_type
) %>%
summarise(Count = n(), .groups = "drop")

carrier_types_plt_permask <- ggplot(part_carriers_summary %>% filter(Carrier_status == "carrier"), aes(x = Mask, y = Count, fill = carrier_type))+ 
    geom_bar(stat="identity", position = "dodge") + facet_wrap(~GENE) + 
    geom_text(stat = "identity", aes(label = Count, color = carrier_type), position = position_dodge(width = 1), vjust = -0.5, size = 3, show.legend = FALSE) +
    theme_minimal() + 
    labs(x = "Mask", y = "Number of Carriers", fill = "Carrier type") + guides(color= guide_legend(label = FALSE)) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position="top")

# Establish the number of non carriers for each gene (static)
# The number of carriers (total), the other graph goes into more detail as to which carrier is which 
carrier_numbers_plt <- ggplot(gene_carriers_all, aes(x= "status", fill=status)) + 
geom_bar(stat="count", position = "dodge") + 
geom_text(stat = "count", aes(label = after_stat(count), color = status), position = position_dodge(width = 0.8), vjust = -0.5, size = 3, show.legend = FALSE) +
theme_minimal() + 
labs(x = "", y = "Count", fill = "") + guides(color= guide_legend(label = FALSE)) + facet_wrap(~GENE)

# Carriers of pathogenic in various different genes (?)
gene_carriers_all %>% 
group_by(SAMPLE) %>% 
summarise(numgenes_carriers=sum(status=='carrier')) %>% pull(numgenes_carriers) %>% table() 

ggsave(filename ="carrier_types_genes_permask_AF.png", carrier_types_plt_permask, width = 12, height = 6, device = "png", dpi = 300)
ggsave(filename= "carrier_numbers_all.png", carrier_numbers_plt,width = 6, height = 6, device = "png", dpi = 300) 

###############################################################################
 
# Transforming (some) variables 

###############################################################################
print(paste("Transforming the sex variable to Female = 0 and Male = 1",
"Transforming the qualifications variable into University = 1 and no University = 0",
"Collapsing the ethnicity categories into White background = 0, Mixed background = 1 and Other = 2",
"Collapsing the genotype arrray to binary array based on batch numbers", collapse = '\n'))

all <- all %>% 
  mutate(sex_coded = ifelse(Sex == "Female", 0, 1), 
         uni_nouni = ifelse(qualifications == 'College or University degree', 1, 0),
         ethnicity_collapsed = case_when(
           ethnicity %in% c('White', 'British','Irish','Any other white background') ~ 0, 
           ethnicity %in% c('Asian or Asian British', 'White and Black African', 'Any other mixed background', 'Mixed' , 'Black or Black British' , 'White and Black Caribbean', 'White and Asian') ~ 1,
           ethnicity %in% c('Chinese', 'Pakistani' , 'African' , 'Do not know' , 'Other ethnic group' , 'Indian' , 'Bangladeshi' , 'Caribbean' , 'Any other black background') ~ 2
         ))
mdd <- mdd %>% filter(MajDepr != -9)
all <- left_join(all, mdd, c("f.eid"= "IID"))

################################################################################

# Summary of the carriers vs non carriers in baseline demographic variables 

################################################################################
# Left join all the carrier status' to the demographic information

all <- reduce(gene_carriers_list, function(x,y){
    left_join(x,y, by = c("f.eid"="SAMPLE"))
}, .init = all)

summary_carrier_status <- function(carrier_dataframe, group_var) {
    carrier_dataframe <- carrier_dataframe %>% filter(!is.na(!!sym(group_var)))
    carrier_summary <- carrier_dataframe %>% 
    group_by(!!sym(group_var)) %>% 
summarise(mean_age = mean(Age, na.rm = T),
sd_age = sd(Age, na.rm = TRUE),
mean_bmi = mean(BMI, na.rm = T),
sd_bmi = sd(BMI, na.rm = TRUE),
mean_TDI = mean(TDI, na.rm = TRUE),
sd_TDI = sd(TDI, na.rm = T),
num_females = sum(sex_coded==0),
num_males = sum(sex_coded == 1),
num_current_smokers = sum(smoking_stat == "Current"),
num_never_smokers = sum(smoking_stat == "Never"),
num_previous_smokers = sum(smoking_stat == "Previous"),
num_uni = sum(uni_nouni == 1),
num_nouni = sum(uni_nouni==0),
num_white_ethnicity = sum(ethnicity_collapsed == 0, na.rm = TRUE),
num_mixed_ethnicity = sum(ethnicity_collapsed == 1, na.rm = TRUE),
num_other_ethnicity = sum(ethnicity_collapsed == 2, na.rm = TRUE),
num_mdd_cases = sum(MajDepr == 2),
num_mdd_controls = sum(MajDepr == 1),
total = n()) %>% 
mutate(Age_stat = paste0(signif(mean_age,3), " (", signif(sd_age, 3), ")"),
BMI_stat = paste0(signif(mean_bmi,3), " (", signif(sd_bmi, 3), ")"),
TDI_stat = paste0(signif(mean_TDI, 3), " (", signif(sd_TDI, 3), ")"),
females_stat = paste0(num_females, " (", signif((num_females/total)*100,3), "%)"),
current_stat = paste0(num_current_smokers, " (", signif((num_current_smokers/total)*100,3), "%)"),
previous_stat = paste0(num_previous_smokers, " (", signif((num_previous_smokers/total)*100,3), "%)"),
never_stat = paste0(num_never_smokers, " (", signif((num_never_smokers/total)*100,3), "%)"),
uni_stat = paste0(num_uni, " (", signif((num_uni/total)*100,3), "%)"),
nouni_stat = paste0(num_nouni, " (", signif((num_nouni/total)*100,3), "%)"),
white_eth_stat = paste0(num_white_ethnicity, " (", signif((num_white_ethnicity/total)*100,3), "%)"),
mixed_eth_stat = paste0(num_mixed_ethnicity, " (", signif((num_mixed_ethnicity/total)*100,3), "%)"),
other_eth_stat = paste0(num_other_ethnicity, " (", signif((num_other_ethnicity/total)*100,3), "%)"),
mdd_cases_stat = paste0(num_mdd_cases, " (", signif((num_mdd_cases/total)*100,3), "%)"),
mdd_controls_stat = paste0(num_mdd_controls, " (", signif((num_mdd_controls/total)*100,3), "%)")
) %>% 
select(!!sym(group_var), total, ends_with("stat")) %>% 
as.data.frame()
return(carrier_summary)
}

genes <- fads_genes$hgnc_symbol
# Use lapply to apply summary_carrier_status() to each gene
lapply(genes, function(gene) {
    # Create the summary for the gene
    carrier_summary <- summary_carrier_status(all, paste0(gene, "_status"))
    
    # Write the summary to a .tsv file
    write.table(carrier_summary, 
                file = paste0(gene, "_allcarriers_summarydemo.tsv"), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
    # Assign to workspace to inspect
    assign(paste0(gene, "_carrier_summary"),carrier_summary)
})

################################################################################

# Chi Squared test of proportions

################################################################################

chi_test <- function(df, gene) {
    table <- table(df[,"MajDepr"], df[, paste0(gene, "_status")])
    rownames(table) <- c("Controls", "Cases")
    chisq <- chisq.test(table)
    return(chisq)
}

FEN1_chires <- chi_test(all, "FEN1")
FADS1_chires <- chi_test(all, "FADS1")
FADS2_chires <- chi_test(all, "FADS2")
FADS3_chires <- chi_test(all, "FADS3")
MYRF_chires <- chi_test(all, "MYRF")
TMEM258_chires <- chi_test(all, "TMEM258")

extract_chi_res <- function(chi_res, gene) {
    statistic <- chi_res$statistic
    df <- chi_res$statistic
    p <- chi_res$p.value
    method <- chi_res$method
    df <- data.frame(gene=gene, statistic=statistic, df=df, p=p, method = method)
    row.names(df) <- ""
    return(df)
}

chi_results_all <- rbind(extract_chi_res(FEN1_chires, "FEN1"), extract_chi_res(FADS1_chires, "FADS1"), extract_chi_res(FADS2_chires, "FADS2"),
extract_chi_res(FADS3_chires, "FADS3"), extract_chi_res(MYRF_chires, "MYRF"), extract_chi_res(TMEM258_chires, "TMEM258"))

chi_res_list <- list(FEN1_chires, FADS1_chires, FADS2_chires, FADS3_chires, MYRF_chires, TMEM258_chires)
extract_chi_values <- function(chi_res, gene){
    observed <- chi_res$observed
    expected <- chi_res$expected
    chi_values_df <- data.frame(gene = gene,
    MDD_status=rownames(observed),
    Carrier_status= c(rep("Carrier",2 ), rep("Non-carrier",2)),
    observed = as.vector(observed),
    expected = as.vector(expected))
    return(chi_values_df)
}
chi_values_list <- mapply(extract_chi_values, chi_res_list, genes, SIMPLIFY = FALSE)
combined_chi_values <- rbindlist(chi_values_list)

write.table(chi_results_all, "all_carriers_mdd_chisq_results.tsv", sep = "\t", row.names =F, quote = F)
write.table(combined_chi_values, "all_carriers_mdd_chisq_exp_obs.tsv", sep = "\t", row.names = F, quote = F)

###############################################################################

# Summary of the carriers vs non carriers in baseline demographic variables 
# In those with metabolomic data at instance 0 

################################################################################

# Read in metabolite cohort info 
metabolite_baseline <- read.csv("data_participant_metabolomics.csv")
# Colnames (some reason are not the same order as the all cohort - be careful 
# "eid" "p31" "p34" "p22189" "p54_i0" "p21000_i0" "p21003_i0" "p21001_i0" "p20116_i0" "p6138_i0"  "p23443_i0" "p23450_i0" "p23444_i0" "p23451_i0" "p23459_i0"
colnames(metabolite_baseline) <- c("f.eid","Sex", "yob", "TDI", "AC", "ethnicity", "Age", "BMI", "smoking_stat", "qualifications",  "f.23443.0.0", "f.23450.0.0", "f.23444.0.0", "f.23451.0.0", "f.23459.0.0")
## Metabolites included are those with significant MR results (DHA the one with MR and colocalisation results)
# f.23443- Degree of Unsaturation
# f.23450 DHA
# f.23444 Omega-3 Fatty Acids
# f.23451 Omega-3 Fatty Acids to TFA
# f.23459 Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio

metabolite_baseline <- metabolite_baseline %>% 
  mutate(sex_coded = ifelse(Sex == "Female", 0, 1), 
         uni_nouni = ifelse(qualifications == 'College or University degree', 1, 0),
         ethnicity_collapsed = case_when(
           ethnicity %in% c('White', 'British','Irish','Any other white background') ~ 0, 
           ethnicity %in% c('Asian or Asian British', 'White and Black African', 'Any other mixed background', 'Mixed' , 'Black or Black British' , 'White and Black Caribbean', 'White and Asian') ~ 1,
           ethnicity %in% c('Chinese', 'Pakistani' , 'African' , 'Do not know' , 'Other ethnic group' , 'Indian' , 'Bangladeshi' , 'Caribbean' , 'Any other black background') ~ 2
         ))


###################################################################################

# Normalise Metabolomic Measures 

###################################################################################

# Metabolite measures names 
nmr_info$matchingmet <- paste('f.', nmr_info$UKB.Field.ID, '.0.0', sep = "")


#### functions to get the metabolite name or the short version of the metabolite name ####

## the full name ##
get_metabolitename <- function(ID_vector) {
  metabolite_names <- c()
  for (i in ID_vector){
    metabolite <- nmr_info$Description[nmr_info$matchingmet == i]
    metabolite_names <- append(metabolite_names, metabolite)
  }
  return(metabolite_names)
}

## the short name ## 
get_metaboliteshort <- function(ID_vector) {
  metabolite_shorts <- c()
  for(i in ID_vector) {
    metabolite <- nmr_info$Biomarker[nmr_info$matchingmet == i]
    metabolite_shorts <- append(metabolite_shorts, metabolite)
  }
  return(metabolite_shorts)
}

# Plot the distributions of the metabolite measures 
meta_measures <- metabolite_baseline %>% 
select(starts_with("f.")) %>% 
rename("Degree of Unsaturation"="f.23443.0.0",
"Docosahexaenoic Acid" = "f.23450.0.0",
"Omega-3 Fatty Acids" ="f.23444.0.0",
"Omega-3 Fatty Acids to Total Fatty Acids percentage"="f.23451.0.0",
"Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio" = "f.23459.0.0"
) %>% 
pivot_longer(., cols = -c("f.eid"), names_to=c("Metabolite"), values_to=c("Value"))

raw_meta_plt <- ggplot(meta_measures, aes(x = Value)) + geom_histogram() +
facet_wrap(~Metabolite) + 
theme_minimal() + 
labs(x = "Raw Metabolite Value", y = "Count", title = paste0("All participants with NMR data: n = ", nrow(metabolite_baseline)))

ggsave(filename ="raw_metabolite_hists.png", raw_meta_plt, width = 10, height = 6, device = "png", dpi = 300)


# Normalise the metabolite values 
norm_metabolite_baseline <- metabolite_baseline %>% 
mutate(f.23443.0.0 = orderNorm(f.23443.0.0)$x.t,
f.23450.0.0 = orderNorm(f.23450.0.0)$x.t,
f.23444.0.0 = orderNorm(f.23444.0.0)$x.t,
f.23451.0.0 = orderNorm(f.23451.0.0)$x.t,
f.23459.0.0 = orderNorm(f.23459.0.0)$x.t)

norm_meta_measures <- norm_metabolite_baseline %>% 
select(starts_with("f.")) %>% 
rename("Degree of Unsaturation"="f.23443.0.0",
"Docosahexaenoic Acid" = "f.23450.0.0",
"Omega-3 Fatty Acids" ="f.23444.0.0",
"Omega-3 Fatty Acids to Total Fatty Acids percentage"="f.23451.0.0",
"Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio" = "f.23459.0.0"
) %>% 
pivot_longer(., cols = -c("f.eid"), names_to=c("Metabolite"), values_to=c("Value"))

norm_meta_plt <- ggplot(norm_meta_measures, aes(x = Value)) + geom_histogram() +
facet_wrap(~Metabolite) + 
theme_minimal() + 
labs(x = "Normalised Metabolite Value", y = "Count", title = paste0("All participants with NMR data: n = ", nrow(metabolite_baseline)))

ggsave(filename ="norm_metabolite_hists.png", norm_meta_plt, width = 10, height = 6, device = "png", dpi = 300)

###################################################################################

# Distribution of metabolomic data per carriers and non carriers 

###################################################################################

norm_metabolite_baseline <- reduce(gene_carriers_list, function(x,y){
    left_join(x,y, by = c("f.eid"="SAMPLE"))
}, .init = metabolite_baseline)
norm_metabolite_baseline <- left_join(norm_metabolite_baseline, mdd, c("f.eid"= "IID"))
summary_metacarrier_status <- function(carrier_dataframe, group_var) {
    carrier_dataframe <- carrier_dataframe %>% filter(!is.na(!!sym(group_var)))
    carrier_summary <- carrier_dataframe %>% 
    group_by(!!sym(group_var)) %>% 
summarise(mean_age = mean(Age, na.rm = T),
sd_age = sd(Age, na.rm = TRUE),
mean_f.23443 = mean(f.23443.0.0, na.rm = TRUE),
sd_f.23443= sd(f.23443.0.0, na.rm = TRUE),
mean_f.23444 = mean(f.23444.0.0, na.rm = TRUE),
sd_f.23444= sd(f.23444.0.0, na.rm = TRUE),
mean_f.23450 = mean(f.23450.0.0, na.rm = TRUE),
sd_f.23450= sd(f.23450.0.0, na.rm = TRUE),
mean_f.23451 = mean(f.23451.0.0, na.rm = TRUE),
sd_f.23451= sd(f.23451.0.0, na.rm = TRUE),
mean_f.23459 = mean(f.23459.0.0, na.rm = TRUE),
sd_f.23459= sd(f.23459.0.0, na.rm = TRUE),
mean_bmi = mean(BMI, na.rm = T),
sd_bmi = sd(BMI, na.rm = TRUE),
mean_TDI = mean(TDI, na.rm = TRUE),
sd_TDI = sd(TDI, na.rm = T),
num_females = sum(sex_coded==0),
num_males = sum(sex_coded == 1),
num_current_smokers = sum(smoking_stat == "Current"),
num_never_smokers = sum(smoking_stat == "Never"),
num_previous_smokers = sum(smoking_stat == "Previous"),
num_uni = sum(uni_nouni == 1),
num_nouni = sum(uni_nouni==0),
num_white_ethnicity = sum(ethnicity_collapsed == 0, na.rm = TRUE),
num_mixed_ethnicity = sum(ethnicity_collapsed == 1, na.rm = TRUE),
num_other_ethnicity = sum(ethnicity_collapsed == 2, na.rm = TRUE),
num_mdd_cases = sum(MajDepr == 2),
num_mdd_controls = sum(MajDepr == 1),
total = n()) %>% 
mutate(Age_stat = paste0(signif(mean_age,3), " (", signif(sd_age, 3), ")"),
BMI_stat = paste0(signif(mean_bmi,3), " (", signif(sd_bmi, 3), ")"),
f.23443.0.0_stat = paste0(signif(mean_f.23443,3), " (", signif(sd_f.23443, 3), ")"),
f.23444.0.0_stat = paste0(signif(mean_f.23444,3), " (", signif(sd_f.23444, 3), ")"),
f.23450.0.0_stat = paste0(signif(mean_f.23450,3), " (", signif(sd_f.23450, 3), ")"),
f.23451.0.0_stat = paste0(signif(mean_f.23451,3), " (", signif(sd_f.23451, 3), ")"),
f.23459.0.0_stat = paste0(signif(mean_f.23459,3), " (", signif(sd_f.23459, 3), ")"),
TDI_stat = paste0(signif(mean_TDI, 3), " (", signif(sd_TDI, 3), ")"),
females_stat = paste0(num_females, " (", signif((num_females/total)*100,3), "%)"),
current_stat = paste0(num_current_smokers, " (", signif((num_current_smokers/total)*100,3), "%)"),
previous_stat = paste0(num_previous_smokers, " (", signif((num_previous_smokers/total)*100,3), "%)"),
never_stat = paste0(num_never_smokers, " (", signif((num_never_smokers/total)*100,3), "%)"),
uni_stat = paste0(num_uni, " (", signif((num_uni/total)*100,3), "%)"),
nouni_stat = paste0(num_nouni, " (", signif((num_nouni/total)*100,3), "%)"),
white_eth_stat = paste0(num_white_ethnicity, " (", signif((num_white_ethnicity/total)*100,3), "%)"),
mixed_eth_stat = paste0(num_mixed_ethnicity, " (", signif((num_mixed_ethnicity/total)*100,3), "%)"),
other_eth_stat = paste0(num_other_ethnicity, " (", signif((num_other_ethnicity/total)*100,3), "%)"),
mdd_cases_stat = paste0(num_mdd_cases, " (", signif((num_mdd_cases/total)*100,3), "%)"),
mdd_controls_stat = paste0(num_mdd_controls, " (", signif((num_mdd_controls/total)*100,3), "%)")
) %>% 
select(!!sym(group_var), total, ends_with("stat")) %>% 
as.data.frame()
return(carrier_summary)
}

# Use lapply to apply summary_metacarrier_status() to each gene
lapply(genes, function(gene) {
    # Create the summary for the gene
    carrier_summary <- summary_metacarrier_status(norm_metabolite_baseline, paste0(gene, "_status"))
    
    # Write the summary to a .tsv file
    write.table(carrier_summary, 
                file = paste0(gene, "_metacarriers_summarydemo.tsv"), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
    # Assign to workspace to inspect
    assign(paste0(gene, "_metacarrier_summary"),carrier_summary)
})

###########################################################################

# Statistically test the differences in metabolite distributions between carriers and non carriers of FADS prioritised variants

###########################################################################
metabolites <- colnames(norm_metabolite_baseline %>% select(ends_with(".0.0")))
run_t_test <- function(metabolite, gene, data) {
    gene_status <- paste0(gene, "_status")
    ttest_result <- t.test(data[[metabolite]]~ data[[gene_status]], data = data)
    result <- data.frame(
        metabolite= get_metabolitename(metabolite),
        gene=gene,
        t_statistic=ttest_result$statistic,
        df=ttest_result$parameter,
        conf_int_lower=ttest_result$conf.int[1],
        conf_int_upper=ttest_result$conf.int[2],
        mean_diff=ttest_result$estimate[1] - ttest_result$estimate[2],
        p_value=ttest_result$p.value
    )
    return(result)
}
ttest_results <- list()
for (metabolite in metabolites) {
    for (gene in genes) {
        ttest_results[[paste(metabolite, gene, sep = "_")]] <- run_t_test(metabolite, gene, norm_metabolite_baseline)
    }
}
ttest_results_df <- do.call(rbind, ttest_results)
ttest_results_df$pval_BH <- p.adjust(ttest_results_df$p_value, method = "BH")
ttest_results_df <- ttest_results_df %>% arrange(pval_BH)

write.table(ttest_results_df, "ttest_metabolite_carrier_all_results.tsv", sep = "\t", row.names = F, quote = F)

###################################################################################

# Plotting the distribution of metabolomic data per carriers and non carriers 

################################################################################### 

data_summary <- function(x) {
  m <- mean(x)
ymin <- m - sd(x)
ymax <- m+sd(x)
  return(c(y=m, ymin = ymin, ymax = ymax))
}

violin_meta_dists <- function(gene, metabolite) {
    # Filter out NA values for the gene status column
    norm_metabolite <- norm_metabolite_baseline %>% 
        filter(!is.na(!!sym(paste0(gene, "_status")))) %>%
        mutate(!!paste0(gene, "_status") := factor(!!sym(paste0(gene, "_status")), levels = c("carrier", "non-carrier")))
    
    # Get the relevant ttest result
    metabolitename <- get_metabolitename(metabolite)
    genename = gene
    ttest_res <- ttest_results_df %>% filter(metabolite == metabolitename & gene == genename)
    ttest_p <- ttest_res$p_value
    ttest_plabel <- paste0("p = ", formatC(ttest_p, format = "e", digits = 2) %>% as.character())
    plt <- ggplot(norm_metabolite, aes(x = as.factor(!!sym(paste0(gene, "_status"))), 
                                       y = !!sym(metabolite),
                                       color = as.factor(!!sym(paste0(gene, "_status"))), 
                                       fill = as.factor(!!sym(paste0(gene, "_status"))))) + 
        geom_violin(trim = FALSE, na.rm = TRUE, alpha = 0.3) +
        scale_color_brewer(palette = "Set1", aesthetics = c("colour", "fill"), 
                           labels = c("Carrier", "Non-carrier"), 
                           name = "Prioritised variant in gene") + 
        theme_classic() +
        geom_jitter(shape = 16, position = position_jitter(0.2), na.rm = TRUE) +
        stat_summary(fun.data = data_summary, shape = 23, color = "black", na.rm = TRUE) + 
        labs(title = "", 
             x = "", 
             y = 'Normalised measure') + 
        scale_x_discrete(labels = c("Carrier", "Non-carrier")) +
        theme(legend.position = "none", 
              text = element_text(size = 9)) +

              geom_signif(comparisons = list(c("carrier", "non-carrier")),
              map_signif_level = TRUE,
              annotations=ttest_plabel,
              tip_length = 0.03, color = "black",
              vjust = 2)

    return(plt)
}

# Plotting the violin plots 
lapply(genes, function(gene) {
      metabolites <- c("f.23443.0.0", "f.23444.0.0", "f.23450.0.0", "f.23451.0.0", "f.23459.0.0")
  
  # Use ggarrange to combine the plots
  plot_list <- lapply(metabolites, function(metabolite) {
    violin_meta_dists(gene, metabolite)
  })
    arranged_plot <- ggarrange(plotlist = plot_list, 
                             nrow = 3, 
                             ncol = 2, 
                             labels = sapply(metabolites, get_metaboliteshort))
  
  # Save the plot
  ggsave(filename = paste0(gene, "_meta_hists_carriers.png"), 
         plot = annotate_figure(arranged_plot, top = text_grob(paste0(gene, ": All prioritised variants"), size = 14, face = "bold")),
         width = 8, height = 10, device = "png", dpi = 300)
})

