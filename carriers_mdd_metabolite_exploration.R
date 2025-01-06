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
library(UpSetR) # Have to install 
library(caret) # Have to install

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
print("Reading in MDD phenotype file")
mdd <- read.table("ukb_unrel_eur_pgc3_mdd.pheno", header = T)
# Filter out missing values for MDD
mdd <- mdd %>% filter(MajDepr != -9)
print("Reading in the Metabolite phenotype file")
metabol <- read.table("ukb_unrel_eur_metabol.pheno", header = T)
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
    mutate(Mask1.singleton_status = ifelse(Mask1.0.01_status == "carrier" & singleton == 1, "carrier", "non-carrier"),
Mask2.singleton_status = ifelse(Mask2.0.01_status == "carrier" & singleton == 1, "carrier", "non-carrier"),
Mask3.singleton_status = ifelse(Mask3.0.01_status == "carrier" & singleton == 1, "carrier", "non-carrier"),
Mask4.singleton_status = ifelse(Mask4.0.01_status == "carrier" & singleton == 1, "carrier", "non-carrier"),
Mask5.singleton_status = ifelse(Mask5.0.01_status == "carrier" & singleton == 1, "carrier", "non-carrier"),
    )
    participant_geno_all <- participant_genotype 
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
part_genotype_all <- rbind(FADS1_part_genotype, FADS2_part_genotype, FADS3_part_genotype, FEN1_part_genotype, MYRF_part_genotype, TMEM258_part_genotype)

# Participant summary information (what type of carriers etc)

part_summary_all_long <- part_summary_all %>% 
pivot_longer(cols = starts_with('Mask'),
names_to="Mask",
values_to = "Carrier_status") %>%
mutate(Mask = gsub("_status", "", Mask))

# Computing separate summaries for the MDD-cohort and the metabolite cohort
part_summary_long_mdd <- part_summary_all_long %>% filter(SAMPLE %in% mdd$IID)
part_summary_long_metabol <- part_summary_all_long %>% filter(SAMPLE %in% metabol$IID)

part_carriers_summary <- part_summary_all_long %>% group_by(
  GENE, Mask, Carrier_status, carrier_type
) %>%
summarise(Count = n(), .groups = "drop")

part_carriers_summary_mdd <- part_summary_long_mdd %>% group_by(
  GENE, Mask, Carrier_status, carrier_type
) %>%
summarise(Count = n(), .groups = "drop")

part_carriers_summary_metabol <- part_summary_long_metabol %>% group_by(
  GENE, Mask, Carrier_status, carrier_type
) %>%
summarise(Count = n(), .groups = "drop")

##### Carrier types per mask for the MDD cohort 
carrier_types_plt_permask_mdd <- ggplot(part_carriers_summary_mdd %>% filter(Carrier_status == "carrier"), aes(x = Mask, y = Count, fill = carrier_type))+ 
    geom_bar(stat="identity", position = "dodge") + facet_wrap(~GENE) + 
    geom_text(stat = "identity", aes(label = Count, color = carrier_type), position = position_dodge(width = 1), vjust = -0.5, size = 3, show.legend = FALSE) +
    theme_minimal() + 
    labs(x = "Mask", y = "Number of Carriers", fill = "Carrier type", title = "MDD-cohort") + guides(color= guide_legend(label = FALSE)) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"))

##### Carrier types per mask for the Metabolite cohort 

carrier_types_plt_permask_metabolite <- ggplot(part_carriers_summary_metabol %>% filter(Carrier_status == "carrier"), aes(x = Mask, y = Count, fill = carrier_type))+ 
    geom_bar(stat="identity", position = "dodge") + facet_wrap(~GENE) + 
    geom_text(stat = "identity", aes(label = Count, color = carrier_type), position = position_dodge(width = 1), vjust = -0.5, size = 3, show.legend = FALSE) +
    theme_minimal() + 
    labs(x = "Mask", y = "Number of Carriers", fill = "Carrier type", title = "Metabolite-cohort") + guides(color= guide_legend(label = FALSE)) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))

#### Saving the plots 
ggsave(filename ="carrier_types_genes_permask_AF_cohorts.png", 
ggarrange(carrier_types_plt_permask_mdd, 
carrier_types_plt_permask_metabolite, 
common.legend = T, nrow = 2, ncol =1, legend = "right"), width = 12, height = 6, device = "png", dpi = 300)

ggsave(filename ="carrier_types_genes_permask_AF_mddcohort.png", 
carrier_types_plt_permask_mdd, width = 12, height = 6, device = "png", dpi = 300)

ggsave(filename ="carrier_types_genes_permask_AF_metabolitecohort.png", 
carrier_types_plt_permask_metabolite, width = 12, height = 6, device = "png", dpi = 300)

#### Plots about gene carriers/non carriers
#### For MDD and Metabolite cohort separately

carrier_summary_all_long <- gene_carriers_all %>% 
pivot_longer(cols = starts_with('Mask'),
names_to="Mask",
values_to = "Carrier_status")  %>%
mutate(Mask = gsub("_status", "", Mask))

gene_carriers_mdd <- gene_carriers_all %>% filter(SAMPLE %in% mdd$IID)
gene_carriers_metabol <- gene_carriers_all %>% filter(SAMPLE %in% metabol$IID)
carrier_summary_all_long_mdd <- carrier_summary_all_long %>% filter(SAMPLE %in% mdd$IID)
carrier_summary_all_long_metabol <- carrier_summary_all_long %>% filter(SAMPLE %in% metabol$IID)

# Establish the number of non carriers for each gene (static)
# The number of carriers (total), the other graph goes into more detail as to which carrier is which 
carrier_numbers_plt <- ggplot(gene_carriers_all, aes(x= "status", fill=status)) + 
geom_bar(stat="count", position = "dodge") + 
geom_text(stat = "count", aes(label = after_stat(count), color = status), position = position_dodge(width = 0.8), vjust = -0.5, size = 2, show.legend = FALSE) +
theme_minimal() + 
labs(x = "", y = "Count", fill = "") + 
guides(color= guide_legend(label = FALSE)) + facet_wrap(~GENE) +
theme(plot.title = element_text(hjust = 0.5, face = "bold"))

carrier_numbers_plt_mdd <- ggplot(gene_carriers_mdd, aes(x= "status", fill=status)) + 
geom_bar(stat="count", position = "dodge") + 
geom_text(stat = "count", aes(label = after_stat(count)),  color = "black",position = position_dodge(width = 0.8), vjust = -0.5, size = 3, show.legend = FALSE) +
theme_minimal() + 
labs(x = "", y = "Count", fill = "", title = "MDD-cohort") + 
guides(color= guide_legend(label = FALSE)) + facet_wrap(~GENE) +
theme(plot.title = element_text(hjust = 0.5, face = "bold"))

carrier_numbers_plt_metabol <- ggplot(gene_carriers_metabol, aes(x= "status", fill=status)) + 
geom_bar(stat="count", position = "dodge") + 
geom_text(stat = "count", aes(label = after_stat(count)), color = "black", position = position_dodge(width = 0.8), vjust = -0.5, size =3, show.legend = FALSE) +
theme_minimal() + 
labs(x = "", y = "Count", fill = "", title = "Metabolite-cohort") + 
guides(color= guide_legend(label = FALSE)) + facet_wrap(~GENE) +
theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(filename ="carrier_numbers_genes_percohort.png", 
ggarrange(carrier_numbers_plt_mdd, carrier_numbers_plt_metabol, common.legend = T, nrow = 2, ncol = 1), 
 width = 21 / 2.54, height = 29.7 / 2.54,device = "png", dpi = 300)

##### Carrier numbers per mask separately for MDD and metabolite cohort 

carrier_numbers_permask <- ggplot(carrier_summary_all_long, aes(x= Mask, fill=Carrier_status)) + 
    geom_bar(stat="count", position = "dodge") + 
    geom_text(stat = "count", aes(label = after_stat(count)), color = "black", position = position_dodge(width = 0.8), vjust = -0.5, size = 2, show.legend = FALSE) +
    theme_minimal() + 
    labs(x = "", y = "Count", fill = "") + guides(color= guide_legend(label = FALSE)) + facet_wrap(~GENE)+theme(
        axis.text.x = element_text(angle = 45, hjust = 1))


carrier_numbers_permask_mdd <- ggplot(carrier_summary_all_long_mdd, aes(x= Mask, fill=Carrier_status)) + 
    geom_bar(stat="count", position = "dodge") + 
    geom_text(stat = "count", aes(label = after_stat(count)), color = "black", position = position_dodge(width = 0.8), vjust = -0.5, size = 2, show.legend = FALSE) +
    theme_minimal() + 
    labs(x = "", y = "Count", fill = "", title= "MDD-cohort") + guides(color= guide_legend(label = FALSE)) + facet_wrap(~GENE)+theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))


carrier_numbers_permask_metabol <- ggplot(carrier_summary_all_long_metabol, aes(x= Mask, fill=Carrier_status)) + 
    geom_bar(stat="count", position = "dodge") + 
    geom_text(stat = "count", aes(label = after_stat(count)), color = "black", position = position_dodge(width = 0.8), vjust = -0.5, size = 2, show.legend = FALSE) +
    theme_minimal() + 
    labs(x = "", y = "Count", fill = "", title = "Metabolite-cohort") + guides(color= guide_legend(label = FALSE)) + facet_wrap(~GENE)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(filename="carrier_numbers_permask_mdd.png", carrier_numbers_permask_mdd, width = 12, height = 8, device = "png", dpi = 300)
ggsave(filename="carrier_numbers_permask_metabolite.png", carrier_numbers_permask_metabol, width = 12, height = 8, device = "png", dpi = 300)

# Carriers of pathogenic in various different genes (?)
# UpSet Plot of MDD cohort cross-carriers 
carriers_across_genes_mdd <- list()
carriers_across_genes_metabol <- list()

gene_carriers_summ_mdd <- gene_carriers_mdd %>% 
group_by(SAMPLE) %>% 
summarise(numgenes_carriers=sum(status=='carrier'),
gene_carriers = paste(GENE[status=="carrier"], collapse = ",")) 

gene_carriers_summ_metabol <- gene_carriers_metabol %>% 
group_by(SAMPLE) %>% 
summarise(numgenes_carriers=sum(status=='carrier'),
gene_carriers = paste(GENE[status=="carrier"], collapse = ",")) 

# List of genes to check
genes <- c("FADS1", "FADS2", "FADS3", "FEN1", "MYRF", "TMEM258")
# Create columns dynamically for each gene
carriers_across_genes_mdd <- genes %>%
  set_names() %>%
  map(~ gene_carriers_summ_mdd %>%
        filter(grepl(.x, gene_carriers)) %>%
        pull(SAMPLE))

carriers_across_genes_metabol <- genes %>%
  set_names() %>%
  map(~ gene_carriers_summ_metabol %>%
        filter(grepl(.x, gene_carriers)) %>%
        pull(SAMPLE))

upset_plt_metabol <- upset(fromList(carriers_across_genes_metabol), nsets = length(carriers_across_genes_metabol), 
set_size.show = TRUE, order.by = 'freq') + grid.text("Metabolite-cohort", x = 0.5, y = 0.97, gp = gpar(fontsize = ))

# Plot the Upset Plot
png(paste0("UpSet_privar_carriers_genes_mdd.png"), 
    width = 3200, height = 2000, res = 300, type = "cairo")
upset(fromList(carriers_across_genes_mdd), nsets = length(carriers_across_genes_mdd),
set_size.show = TRUE, order.by = 'freq')
grid.text("MDD-cohort", x = 0.5, y = 0.97, gp = gpar(fontsize = 20))
dev.off()

png(paste0("UpSet_privar_carriers_genes_metabol.png"), 
    width = 3200, height = 2000, res = 300, type = "cairo")
upset(fromList(carriers_across_genes_metabol), nsets = length(carriers_across_genes_metabol),
set_size.show = TRUE, order.by = 'freq')
grid.text("Metabolite-cohort", x = 0.5, y = 0.97, gp = gpar(fontsize = 20))
dev.off()

ggsave(filename ="carrier_types_genes_permask_AF.png", carrier_types_plt_permask, width = 12, height = 6, device = "png", dpi = 300)
ggsave(filename= "carrier_numbers_all.png", carrier_numbers_plt,width = 6, height = 8, device = "png", dpi = 300) 
ggsave(filename="carrier_numbers_permask.png", carrier_numbers_permask, width = 10, height = 8, device = "png", dpi = 300)

###############################################################################
 
# Transforming (some) variables 

###############################################################################

print("The covariate file")
covars <- read.table("ukb_unrel_eur_covars.covar", header = T)
# Read in all cohort info 
all <- read.csv("data_participant_all.csv")
colnames(all) <- c("f.eid", "Sex", "yob", "TDI", "AC", "ethnicity", "Age","BMI", "smoking_stat", "qualifications",
 "f.23444.0.0", "f.23451.0.0", "f.23459.0.0", "f.23443.0.0", "f.23450.0.0",
 "nmr_proc_batch","spectrometer", "nmr_shipment_plate",
 paste0("PC", 1:20), "genotyping_array_batch"
 )
print(paste("Transforming the sex variable to Female = 0 and Male = 1",
"Transforming the qualifications variable into University = 1 and no University = 0",
"Collapsing the ethnicity categories into White background = 0, Mixed background = 1 and Other = 2",
"Collapsing the genotype arrray to binary array based on batch numbers",
"Transforming the genotype array variable to BiLEVE = 1 and Axiom = 0", collapse = '\n'))

all <- all %>% 
  mutate(sex_coded = ifelse(Sex == "Female", 0, 1), 
         uni_nouni = ifelse(qualifications == 'College or University degree', 1, 0),
         genotype_array = case_when(grepl("BiLEVE", genotyping_array_batch) ~ 1, 
         grepl("Batch", genotyping_array_batch) ~ 0,
         TRUE ~ NA_real_),
         ethnicity_collapsed = case_when(
           ethnicity %in% c('White', 'British','Irish','Any other white background') ~ 0, 
           ethnicity %in% c('Asian or Asian British', 'White and Black African', 'Any other mixed background', 'Mixed' , 'Black or Black British' , 'White and Black Caribbean', 'White and Asian') ~ 1,
           ethnicity %in% c('Chinese', 'Pakistani' , 'African' , 'Do not know' , 'Other ethnic group' , 'Indian' , 'Bangladeshi' , 'Caribbean' , 'Any other black background') ~ 2
         ))

all <- left_join(all, mdd, c("f.eid"= "IID"))

################################################################################

# Summary of the carriers vs non carriers in baseline demographic variables 

################################################################################
# Left join all the carrier status' to the demographic information

# Mutate the dataframe to classify non-carriers across all the columns (not just a particular mask)
# Define the mask columns
mask_columns <- grep("Mask", names(gene_carriers_all), value = TRUE)

# Apply the transformation
gene_carriers_mdd <- gene_carriers_mdd %>%
  mutate(across(
    all_of(mask_columns), 
    ~ ifelse(status == "non-carrier", "non_carrier_any", .)
  ))
gene_carriers_metabol <- gene_carriers_metabol %>%
  mutate(across(
    all_of(mask_columns), 
    ~ ifelse(status == "non-carrier", "non_carrier_any", .)
  ))


# Merge the general carrier information with the demographics (i.e carrier/non_carrier)
all_gene_carriers_mdd <- left_join(all, gene_carriers_mdd, by = c("f.eid"="SAMPLE"))
all_gene_carriers_metabol <- left_join(all, gene_carriers_metabol, by = c("f.eid"="SAMPLE"))

summary_carrier_status <- function(carrier_dataframe, gene) {
  carrier_dataframe <- carrier_dataframe %>% filter(GENE == gene)
    carrier_dataframe <- carrier_dataframe %>% filter(!is.na(status))
    carrier_summary <- carrier_dataframe %>% 
    group_by(status) %>% 
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
select(status, total, ends_with("stat")) %>% 
as.data.frame()
return(carrier_summary)
}

genes <- fads_genes$hgnc_symbol
# Use lapply to apply summary_carrier_status() to each gene
lapply(genes, function(gene) {
    # Create the summary for the gene
    print(gene)
    mddcohort_carrier_summary <- summary_carrier_status(all_gene_carriers_mdd, gene)
    metabolcohort_carrier_summary <- summary_carrier_status(all_gene_carriers_metabol, gene)
    # Write the summary to a .tsv file
    write.table(mddcohort_carrier_summary, 
                file = paste0(gene, "_allcarriers_summarydemo_MDDcohort.tsv"), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
    write.table(metabolcohort_carrier_summary, 
                file = paste0(gene, "_allcarriers_summarydemo_metabolitecohort.tsv"), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
    # Assign to workspace to inspect
    assign(paste0(gene, "_carrier_summary_mdd"),mddcohort_carrier_summary)
    assign(paste0(gene, "_carrier_summary_metabol"),metabolcohort_carrier_summary)
})

################################################################################

# Chi Squared test of proportions (now not needed as have run REGENIE instead)

################################################################################
all_gene_carriers_mdd <- all_gene_carriers_mdd %>%
  mutate(across(
    all_of(mask_columns), 
    ~ ifelse(. == "non_carrier_any", "non-carrier", .)
  ))

# Extract Chi-square results
extract_chi_res <- function(chi_res, gene, mask) {
    statistic <- chi_res$statistic
    df <- chi_res$parameter
    p <- chi_res$p.value
    method <- chi_res$method
    df <- data.frame(gene=gene, mask=mask, statistic=statistic, df=df, p=p, method = method)
    row.names(df) <- ""
    return(df)
}
# Extract observed and expected values from Chi-square test
extract_chi_values <- function(chi_res, gene, mask){
    observed <- chi_res$observed
    expected <- chi_res$expected
    chi_values_df <- data.frame(gene = rep(gene,4),
    mask = rep(mask, 4),
    MDD_status=rep(rownames(observed),2),
    Carrier_status= c(rep("Carrier",2 ), rep("Non-carrier",2)),
    observed = as.vector(observed),
    expected = as.vector(expected))
    return(chi_values_df)
}

# Main chi-square test function
chi_test <- function(df, gene) {
  # Filter data for the given gene
  df <- df %>% filter(GENE == gene)
  
  # Define the masks to process
  masks <- c("Mask2.0.01", "Mask3.0.01", "Mask4.0.01", "Mask5.0.01")
  singleton_masks <- c("Mask2.singleton", "Mask3.singleton", "Mask4.singleton", "Mask5.singleton")
  
  # Include Mask1 only if gene is not FEN1
  if (gene != "FEN1") {
    masks <- c("Mask1.0.01", masks)
    singleton_masks <- c("Mask1.singleton", singleton_masks)
  }
  
  # Helper function to create tables, run chi-sq tests, and extract results
  process_masks <- function(mask_list, prefix) {
    results <- lapply(mask_list, function(mask) {
      table_data <- table(df[,"MajDepr"], df[[paste0(mask, "_status")]])
      rownames(table_data) <- c("Controls", "Cases")
      chi_test <- chisq.test(table_data)
      
      list(
        res = extract_chi_res(chi_test, gene, paste0(prefix, mask)),
        values = extract_chi_values(chi_test, gene, paste0(prefix, mask))
      )
    })
    results
  }
  
  # Process both sets of masks
  chi_mask_results <- process_masks(masks, "")
  chi_singleton_results <- process_masks(singleton_masks, "")
  
  # Combine all results
  chi_results_all_masks <- do.call(rbind, lapply(chi_mask_results, `[[`, "res"))
  chi_results_all_masks <- rbind(chi_results_all_masks, do.call(rbind, lapply(chi_singleton_results, `[[`, "res")))
  
  chi_res_values_all <- do.call(rbind, lapply(chi_mask_results, `[[`, "values"))
  chi_res_values_all <- rbind(chi_res_values_all, do.call(rbind, lapply(chi_singleton_results, `[[`, "values")))
  
  return(list(chi_results_all_masks, chi_res_values_all))
}


FEN1_chires <- chi_test(all_gene_carriers_mdd, "FEN1")
FADS1_chires <- chi_test(all_gene_carriers_mdd, "FADS1")
FADS2_chires <- chi_test(all_gene_carriers_mdd, "FADS2")
FADS3_chires <- chi_test(all_gene_carriers_mdd, "FADS3")
MYRF_chires <- chi_test(all_gene_carriers_mdd, "MYRF")
TMEM258_chires <- chi_test(all_gene_carriers_mdd, "TMEM258")

chi_results_all <- rbind(FADS1_chires[[1]] %>% as.data.frame(), FADS2_chires[[1]] %>% as.data.frame(), FADS3_chires[[1]] %>% as.data.frame(), FEN1_chires[[1]] %>% as.data.frame(), MYRF_chires[[1]] %>% as.data.frame(), TMEM258_chires[[1]] %>% as.data.frame()) 
chi_values_all <- rbind(FADS1_chires[[2]], FADS2_chires[[2]], FADS3_chires[[2]], FEN1_chires[[2]], MYRF_chires[[2]], TMEM258_chires[[2]]) 

write.table(chi_results_all, "all_carriers_mdd_chisq_results_permask.tsv", sep = "\t", row.names =F, quote = F)
write.table(chi_values_all, "all_carriers_mdd_chisq_exp_obs_permask.tsv", sep = "\t", row.names = F, quote = F)

###############################################################################

# Confusion matrices for each Mask and MDD 

###############################################################################

# Plot the confusion matrices between the carriers of a particular mask and non-carriers of that particular mask (these non carriers could be carriers of prioritised variants in another mask)

chi_values_long <- reshape2::melt(chi_values_all, id.vars = c("gene", "mask", "MDD_status", "Carrier_status"), 
                            measure.vars = c("observed", "expected"), 
                            variable.name = "Type", value.name = "Value")
chi_values_long <- chi_values_long %>% mutate(Type = case_when(Type == "observed" ~ "Observed",
Type == "expected"~ "Expected"), fill=paste0(MDD_status, ":", Carrier_status))


chi_exp_obs_plt <- function(genename) {
plt <- ggplot(chi_values_long %>% filter(gene == genename), aes(x = Type, y = MDD_status, fill = Value)) +
geom_tile(color = "black") +
geom_text(aes(label = round(Value, 1)), color = "black") +
facet_grid(mask ~ Carrier_status) +
scale_fill_gradient(name = "Count",low = "white", high = "red") +
labs(
title = paste0(genename, ": Observed and expected MDD status and carrier status \ncounts per mask"),
x = "Carrier Status",
y = "MDD Status",
fill = "Value"
) +
theme_minimal()+
 theme(
strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5, face = "bold"))
return(plt)
}

ggsave(filename ="chi_exp_obs_FADS1.png", chi_exp_obs_plt("FADS1"), width = 8, height = 8, device = "png", dpi = 300)
ggsave(filename ="chi_exp_obs_FADS2.png", chi_exp_obs_plt("FADS2"), width = 8, height = 8, device = "png", dpi = 300)
ggsave(filename ="chi_exp_obs_FADS3.png", chi_exp_obs_plt("FADS3"), width = 8, height = 8, device = "png", dpi = 300)
ggsave(filename ="chi_exp_obs_FEN1.png", chi_exp_obs_plt("FEN1"), width = 8, height = 8, device = "png", dpi = 300)
ggsave(filename ="chi_exp_obs_MYRF.png", chi_exp_obs_plt("MYRF"), width = 8, height = 8, device = "png", dpi = 300)
ggsave(filename ="chi_exp_obs_TMEM258.png", chi_exp_obs_plt("TMEM258"), width = 8, height = 8, device = "png", dpi = 300)

observed_vals_plt <- ggplot(chi_values_long %>% filter(Type == "Observed"), aes(x = Carrier_status, y = MDD_status, fill = Value)) +
geom_tile(color = "black") +
geom_text(aes(label = round(Value, 1)), color = "black", size = 2.5) +
facet_grid(mask~gene) +
scale_fill_gradient(name = "Count",low = "white", high = "red") +
labs(
title = paste0("Observed counts of MDD and carrier status"),
x = "Carrier Status",
y = "MDD Status",
fill = "Value"
) +
 theme_minimal()+
theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 45, hjust = 1),
plot.title = element_text(hjust = 0.5, face = "bold"))

expected_vals_plt <- ggplot(chi_values_long %>% filter(Type == "Expected"), aes(x = Carrier_status, y = MDD_status, fill = Value)) +
geom_tile(color = "black") +
geom_text(aes(label = round(Value, 0)), color = "black", size = 2.5) +
facet_grid(mask~gene) +
scale_fill_gradient(name = "Count",low = "white", high = "red") +
labs(
title = paste0("Expected counts of MDD and carrier status"),
x = "Carrier Status",
y = "MDD Status",
fill = "Value"
) +
 theme_minimal()+
theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 45, hjust = 1),
plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(filename ="chi_observed_allgenes.png", observed_vals_plt, width = 8, height = 8, device = "png", dpi = 300)
ggsave(filename ="chi_expected_allgenes.png", expected_vals_plt, width = 8, height = 8, device = "png", dpi = 300)


# Summary of the carriers vs non carriers in baseline demographic variables 
# In those with metabolomic data at instance 0 

################################################################################

metabolite_baseline <- all %>% filter(spectrometer != "")

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
distinct(f.eid, .keep_all=TRUE) %>%
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
labs(x = "Raw Metabolite Value", y = "Count", title = paste0("All participants with NMR data: n = ", nrow(metabolite_baseline %>% 
distinct(f.eid, .keep_all=TRUE) )))

ggsave(filename ="raw_metabolite_hists.png", raw_meta_plt, width = 10, height = 6, device = "png", dpi = 300)


# Normalise the metabolite values 
norm_metabolite_baseline <- metabolite_baseline %>% 
mutate(f.23443.0.0 = orderNorm(f.23443.0.0)$x.t,
f.23450.0.0 = orderNorm(f.23450.0.0)$x.t,
f.23444.0.0 = orderNorm(f.23444.0.0)$x.t,
f.23451.0.0 = orderNorm(f.23451.0.0)$x.t,
f.23459.0.0 = orderNorm(f.23459.0.0)$x.t)

norm_meta_measures <- norm_metabolite_baseline %>% 
distinct(f.eid, .keep_all=TRUE) %>%
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
labs(x = "Normalised Metabolite Value", y = "Count", title = paste0("All participants with NMR data: n = ", nrow(metabolite_baseline %>% 
distinct(f.eid, .keep_all=TRUE))))

ggsave(filename ="norm_metabolite_hists.png", norm_meta_plt, width = 10, height = 6, device = "png", dpi = 300)

###################################################################################

# Distribution of metabolomic data per carriers and non carriers 

###################################################################################

norm_metabolite_carriers <- left_join(norm_metabolite_baseline, gene_carriers_all, by = c("f.eid"="SAMPLE"))

summary_metacarrier_status <- function(carrier_dataframe, gene) {
  print(gene)
  carrier_dataframe <- carrier_dataframe %>% filter(GENE == gene)
    carrier_dataframe <- carrier_dataframe %>% filter(!is.na(status))
    carrier_summary <- carrier_dataframe %>% 
    group_by(status) %>% 
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
select(status, total, ends_with("stat")) %>% 
as.data.frame()
return(carrier_summary)
}

# Use lapply to apply summary_metacarrier_status() to each gene
lapply(genes, function(gene) {
    # Create the summary for the gene
    carrier_summary <- summary_metacarrier_status(norm_metabolite_carriers, gene)
    
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

metabolites <- colnames(norm_metabolite_carriers %>% select(ends_with(".0.0")))
norm_metabolite_carriers <- norm_metabolite_carriers %>%
  mutate(across(
    all_of(mask_columns), 
    ~ ifelse(. == "non_carrier_any", "non-carrier", .)
  ))

run_t_test <- function(metabolite, gene, mask, data) {
  print(gene)
  print(mask)
  print(metabolite)
    data <- data %>% filter(GENE == gene)

    if(length(table(data[[mask]]))==1) {
      message(paste0("Skipping:", gene, mask, " - Mask not available"))
      return(NULL)
    }

    print("Running t-test")
    ttest_result <- t.test(data[[metabolite]]~ data[[mask]], data = data)
    print("Results table")
    result <- data.frame(
        metabolite= get_metabolitename(metabolite),
        gene=gene,
        mask=mask,
        t_statistic=ttest_result$statistic,
        df=ttest_result$parameter,
        conf_int_lower=ttest_result$conf.int[1],
        conf_int_upper=ttest_result$conf.int[2],
        mean_diff=ttest_result$estimate[1] - ttest_result$estimate[2],
        p_value=ttest_result$p.value
    )
    return(result)
}

ttest_results <- do.call(
  rbind,
  lapply(metabolites, function(metabolite) {
    do.call(
      rbind, 
      lapply(genes, function(gene) {
        do.call(rbind,
        lapply(mask_columns, function(mask) {
          run_t_test(metabolite, gene, mask, norm_metabolite_carriers)
        })
        )
      })
    )
  })
)

ttest_results <- ttest_results[!sapply(ttest_results, is.null),]

ttest_results$pval_BH <- p.adjust(ttest_results$p_value, method = "BH")
ttest_results <- ttest_results %>% arrange(pval_BH)

write.table(ttest_results, "ttest_metabolite_carrier_permask_all_results.tsv", sep = "\t", row.names = F, quote = F)

###################################################################################

# Plotting the distribution of metabolomic data per carriers and non carriers 

################################################################################### 

data_summary <- function(x) {
  m <- mean(x)
ymin <- m - sd(x)
ymax <- m+sd(x)
  return(c(y=m, ymin = ymin, ymax = ymax))
}

violin_meta_dists <- function(gene, metabolite, mask) {
    # Filter out NA values for the gene status column
    norm_metabolite <- norm_metabolite_carriers %>% 
        filter(GENE == gene) %>%
        mutate(across(
    all_of(mask_columns), 
    ~ factor(., levels = c("carrier", "non-carrier"))
  ))
    # Get the relevant ttest result
    metabolitename <- get_metabolitename(metabolite)
    genename <- gene
    maskname <- mask
    ttest_res <- ttest_results %>% filter(metabolite == metabolitename & gene == genename & mask == maskname)
    ttest_p <- ttest_res$p_value
    ttest_plabel <- paste0("p = ", formatC(ttest_p, format = "e", digits = 2) %>% as.character())
    plt <- ggplot(norm_metabolite, aes(x = as.factor(!!sym(mask)), 
                                       y = !!sym(metabolite),
                                       color = as.factor(!!sym(mask)), 
                                       fill = as.factor(!!sym(mask)))) + 
        geom_violin(trim = FALSE, na.rm = TRUE, alpha = 0.3) +
        scale_color_brewer(palette = "Set1", aesthetics = c("colour", "fill"), 
                           labels = c("Carrier", "Non-carrier"), 
                           name = "Prioritised variant in gene") + 
        theme_classic() +
        geom_jitter(shape = 16, position = position_jitter(0.2), na.rm = TRUE) +
        stat_summary(fun.data = data_summary, shape = 23, color = "black", na.rm = TRUE) + 
        labs(title ="", 
             x = mask, 
             y =  'Normalised measure') + 
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
  lapply(mask_columns, function(mask) {
plot_list <- lapply(metabolites, function(metabolite) {
    violin_meta_dists(gene, metabolite, mask)
  })
  # Use ggarrange to combine the plots
    arranged_plot <- ggarrange(plotlist = plot_list, 
                             nrow = 3, 
                             ncol = 2, 
                             labels = sapply(metabolites, get_metaboliteshort))
  
  # Save the plot
  ggsave(filename = paste0(gene, "_meta_hists_carriers_", gsub("_status", "", mask), ".png"), 
         plot = annotate_figure(arranged_plot, top = text_grob(paste0(gene, ":", gsub("_status", "", mask)), size = 14, face = "bold")),
         width = 8, height = 10, device = "png", dpi = 300)
})
  })

########## Getting variant summaries for the prioritised variants in the MDD cohort and metabolite cohort

# Read in the priority variant tables 
# Subset to the variants which are present in the participants in the two different cohorts
# Save for making the variant plots separately for the MDD cohort and the metabolite cohort 

for (i in c(1:nrow(fads_genes))){
  gene <- fads_genes$hgnc_symbol[i]
  priority <- read.table(paste0(gene, "_priority_annot.tsv"), sep = "\t", header = T)
  print(paste0("Read in priority variants: ", gene, " (n=", nrow(priority), ")"))
  part_genotype_mdd <- part_genotype_all %>% filter(SAMPLE %in% mdd$IID)
  part_genotype_metabol <- part_genotype_all %>% filter(SAMPLE %in% metabol$IID)
  mdd_variants <- c(part_genotype_mdd$het_variants, part_genotype_mdd$hom_variants) %>% unique()
  metabol_variants <- c(part_genotype_metabol$het_variants, part_genotype_metabol$hom_variants) %>% unique()
  priority_mdd <- priority %>% filter(chr_pos_ref_alt %in% mdd_variants)
  priority_metabol <- priority %>% filter(chr_pos_ref_alt %in% metabol_variants)
  write.table(priority_mdd, paste0(gene, "_priority_var_mdd.tsv"), sep = "\t", row.names = F, quote = F)
  write.table(priority_metabol, paste0(gene, "_priority_var_metabol.tsv"), sep = "\t", row.names = F, quote = F)
  assign(paste0(gene, "_priority_var_mdd"), priority_mdd)
  assign(paste0(gene, "_priority_var_metabol"), priority_metabol)
}

