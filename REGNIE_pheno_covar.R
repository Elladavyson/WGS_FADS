
library(data.table)
library(dplyr)
library(bestNormalize) # Have to install
library(ukbnmr) # Have to install
library(purrr)
library(rlang)
library(tidyr)
library(ggplot2)

# dx download /Input/MajorDepression.ukb24262.2021-07.txt
# dx download /Input/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id
# dx download data_participant_all.csv
# dx download data_participant_metabolomics.csv
# dx download /Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3087_v1_norm.vcf.gz

print("Reading in MDD phenotype")
mdd <- read.table("MajorDepression.ukb24262.2021-07.txt", header = T)
print("Reading in cohort data")
# Read in all cohort info 
all <- read.csv("data_participant_all.csv")
# colnames demographics : "eid" "p31"  "p34" "p22189" "p54_i0" "p21000_i0" "p21003_i0" "p21001_i0" "p20116_i0" "p6138_i0"  
# metabolites : "p23444_i0" "p23451_i0" "p23459_i0" "p23443_i0" "p23450_i0"
# f.23443- Degree of Unsaturation
# f.23450 DHA
# f.23444 Omega-3 Fatty Acids
# f.23451 Omega-3 Fatty Acids to TFA
# f.23459 Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio
# metabolite preprocessing variables: p20282_i0, p23650_i0, p23649_i0
# Genotyping covariates: p22009_a{1-20}, p22000
colnames(all) <- c("f.eid", "Sex", "yob", "TDI", "AC", "ethnicity", "Age","BMI", "smoking_stat", "qualifications",
 "f.23444.0.0", "f.23451.0.0", "f.23459.0.0", "f.23443.0.0", "f.23450.0.0",
 "nmr_proc_batch","spectrometer", "nmr_shipment_plate",
 paste0("PC", 1:20), "genotyping_array_batch"
 )
print(paste("Transforming the sex variable to Female = 0 and Male = 1",
"Transforming the qualifications variable into University = 1 and no University = 0",
"Collapsing the ethnicity categories into White background = 0, Mixed background = 1 and Other = 2",
"Collapsing the genotype arrray to binary array based on batch numbers", collapse = '\n'))


all <- all %>% 
  mutate(sex_coded = ifelse(Sex == "Female", 0, 1), 
         uni_nouni = ifelse(qualifications == 'College or University degree', 1, 0),
         genotype_array = case_when(grepl("BiLEVE", genotyping_array_batch) ~ "BiLEVE", 
         grepl("Batch", genotyping_array_batch) ~ "Axiom",
         TRUE ~ "NA"),
         ethnicity_collapsed = case_when(
           ethnicity %in% c('White', 'British','Irish','Any other white background') ~ 0, 
           ethnicity %in% c('Asian or Asian British', 'White and Black African', 'Any other mixed background', 'Mixed' , 'Black or Black British' , 'White and Black Caribbean', 'White and Asian') ~ 1,
           ethnicity %in% c('Chinese', 'Pakistani' , 'African' , 'Do not know' , 'Other ethnic group' , 'Indian' , 'Bangladeshi' , 'Caribbean' , 'Any other black background') ~ 2
         ))
all <- left_join(all, mdd, c("f.eid"= "IID"))
# Replacing blank values as NA for REGENIE
all[all==""]<- NA
lapply(all %>% select(starts_with("f.")), function(column) table(is.na(column)))

###############################################

# Normalise the metabolite data 

###############################################

all <- all %>% 
mutate(f.23443.0.0 = orderNorm(f.23443.0.0)$x.t,
f.23450.0.0 = orderNorm(f.23450.0.0)$x.t,
f.23444.0.0 = orderNorm(f.23444.0.0)$x.t,
f.23451.0.0 = orderNorm(f.23451.0.0)$x.t,
f.23459.0.0 = orderNorm(f.23459.0.0)$x.t)

norm_meta_measures <- all %>% 
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
labs(x = "Normalised Metabolite Value", y = "Count", title = paste0("All participants with NMR data: n=", nrow(all %>% filter(!is.na(spectrometer)))))


# Read in the unrelated participants of european ancestry
unrelated_eur <- read.table("ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id", header = F)
colnames(unrelated_eur) <- c("IID", "FID")
all_unrelated <- all %>% filter(f.eid %in% unrelated_eur$IID)

# Make sure they have WGS data
# bcftools query -l ukb24310_c11_b3087_v1_norm.vcf.gz > sample_names.txt
wgs_samples <- readr::read_lines('sample_names.txt')
all_unrelated_wgs <- all_unrelated %>% 
filter(f.eid %in% wgs_samples)

## Distribution of metabolite levels in this sample (still normal)

norm_meta_measures_unrel_eur <- all_unrelated_wgs %>% 
select(starts_with("f.")) %>% 
rename("Degree of Unsaturation"="f.23443.0.0",
"Docosahexaenoic Acid" = "f.23450.0.0",
"Omega-3 Fatty Acids" ="f.23444.0.0",
"Omega-3 Fatty Acids to Total Fatty Acids percentage"="f.23451.0.0",
"Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio" = "f.23459.0.0"
) %>% 
pivot_longer(., cols = -c("f.eid"), names_to=c("Metabolite"), values_to=c("Value"))

norm_meta_plt_unrel_eur <- ggplot(norm_meta_measures_unrel_eur, aes(x = Value)) + geom_histogram() +
facet_wrap(~Metabolite) + 
theme_minimal() + 
labs(x = "Normalised Metabolite Value", y = "Count", title = paste0("All unrelated participants of European ancestry with NMR data: n=", nrow(all_unrelated_wgs %>% filter(!is.na(spectrometer)))))


######################

# Phenotype file

# All the metabolites and MDD

######################

mdd_pheno <- all_unrelated_wgs %>% 
select("f.eid", "FID", "MajDepr") %>%
rename("IID"="f.eid") %>%
select(FID, IID, MajDepr)
# replace all missing values (coded as -9) with NA 
mdd_pheno$MajDepr[mdd_pheno$MajDepr == -9] <- NA
lapply(mdd_pheno %>% select(-c(FID, IID)), function(column) table(column, useNA="always"))

metabol_pheno <- all_unrelated_wgs %>% 
filter(!is.na(spectrometer)) %>% 
select(starts_with("f."), "FID") %>% 
rename("IID"="f.eid") %>%
select(FID, IID, starts_with("f."))
lapply(metabol_pheno %>% select(starts_with("f.")), function(column) table(is.na(column)))

write.table(mdd_pheno, "ukb_unrel_eur_pgc3_mdd.pheno", sep = "\t", row.names = F, quote = F)
write.table(metabol_pheno, "ukb_unrel_eur_metabol.pheno", sep = "\t", row.names = F, quote=F)

########################

# Covariate files 

########################

covars <- all_unrelated_wgs %>% 
select(f.eid, FID, Age, sex_coded, genotype_array, AC, spectrometer, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
rename("IID"="f.eid") %>% 
select(FID, IID, everything())
lapply(covars %>% select(-c(FID, IID, starts_with("PC"))), function(column) table(column, useNA="always"))

write.table(covars, "ukb_unrel_eur_covars.covar", sep = "\t", row.names = F, quote=F)

# dx upload ukb_unrel_eur_covars.covar
# dx upload ukb_unrel_eur_pgc3_mdd.pheno
# dx upload ukb_unrel_eur_metabol.pheno