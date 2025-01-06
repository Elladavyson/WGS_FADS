#### Read in the demographic table 

mdd_demo <- read.table("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Intermediate_data_files/mdd3_pgc_eur_unrel_sample_demo.tsv", sep = "\t", header = T)
metabol_demo <- read.table("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Intermediate_data_files/baseline_metabolite_eur_unrel_sample_demo.tsv", sep = "\t", header = T)

mdd_demo <- mdd_demo %>% mutate(cohort = "MDD cohort")
metabol_demo <- metabol_demo %>% mutate(cohort = "Metabolite cohort")

process_demo <- function(table) {
  table <- table %>% mutate(smoking = paste0(current_stat, ":", previous_stat, ":", never_stat),
                            ethnicity = paste0(white_eth_stat, ":", mixed_eth_stat, ":", other_eth_stat),
                            mdd = paste0(mdd_cases_stat, ":", mdd_controls_stat)) %>% rename(
                              'Age' = 'Age_stat',
                              'BMI' = 'BMI_stat',
                              'SES' = "TDI_stat", 
                              '% Female' = "females_stat",
                              '% Never smoked' = 'never_stat',
                              '% Attended University/College' = 'uni_stat',
                              "Smoking status" = 'smoking',
                              'MDD' = 'mdd', 
                              'Cohort' = 'cohort',
                              'Total' = 'total',
                              'Ethnicity'= 'ethnicity'
                            )
  table <- table %>% 
    select(Cohort, Total, Age, BMI, SES, `% Female`, `% Attended University/College`, `Smoking status`, Ethnicity, MDD)
  return(table)
}

demo_both <- rbind(mdd_demo, metabol_demo)
demo_both <- demo_both %>% select(cohort, everything())
View(process_demo(demo_both))

