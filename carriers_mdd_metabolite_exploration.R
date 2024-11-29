library(data.table)
library(dplyr)
library(bestNormalize)

# carriers info 
print("Reading in carrier information")
FEN1_carriers <- read.table("FEN1_carrier_info.tsv", sep = "\t", header = T)
mdd <- read.table("MajorDepression.ukb24262.2021-07.txt", header = T)

# Read in all cohort info 
all <- read.csv("data_participant_all.csv")
# colnames "eid" "p31"  "p34" "p22189" "p54_i0" "p21000_i0" "p21003_i0" "p21001_i0" "p20116_i0" "p6138_i0"  "p23449_i0" "p23444_i0" "p23450_i0" "p23457_i0" "p23459_i0" "p23443_i0"
colnames(all) <- c("f.eid", "Sex", "yob", "TDI", "AC", "ethnicity", "Age","BMI", "smoking_stat", "qualifications", "f.23449.0.0", "f.23444.0.0", "f.23450.0.0", "f.23457.0.0", "f.23459.0.0", "f.23443.0.0")

extra_mets <- all %>% filter(f.eid %in% metabolite_baseline$f.eid) %>% select(f.23459.0.0, f.23443.0.0)

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

all <- left_join(all, mdd, c("f.eid"= "IID"))

################################################################################

# Summary of the carriers vs non carriers in baseline demographic variables 

################################################################################

all <- left_join(all, FEN1_carriers, by=c("f.eid"="SAMPLE"))

# Filter out those without FEN1 variant carrier status 
all_FEN1 <- all %>% filter(!is.na(status))

summary_carrier_status <- function(carrier_dataframe) {
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

FEN1_carrier_summary <- summary_carrier_status(all_FEN1)

################################################################################

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

metabolite_baseline <- left_join(metabolite_baseline, mdd, c("f.eid"= "IID"))
metabolite_baseline <- left_join(metabolite_baseline, FEN1_carriers, by=c("f.eid"="SAMPLE"))
meta_bl_FEN1 <- metabolite_baseline %>% filter(!is.na(status))

FEN1_metabolite_carrier_summary <- summary_carrier_status(meta_bl_FEN1)

###################################################################################

# Normalise Metabolomic Measures 

###################################################################################

