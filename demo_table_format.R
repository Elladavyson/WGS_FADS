#### Read in the demographic table 
fads_genes <- read.table("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Intermediate_data_files/FADS_cluster_UKB_pVCF.tsv", sep= "\t", header=T)
mdd_demo <- read.table("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Intermediate_data_files/mdd3_pgc_eur_unrel_sample_demo.tsv", sep = "\t", header = T)
metabol_demo <- read.table("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Intermediate_data_files/baseline_metabolite_eur_unrel_sample_demo.tsv", sep = "\t", header = T)

mdd_demo <- mdd_demo %>% mutate(cohort = "MDD cohort")
metabol_demo <- metabol_demo %>% mutate(cohort = "Metabolite cohort")

process_demo <- function(table, carrier_status) {
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
  if(carrier_status == TRUE) {
    table <- table %>% 
      mutate(status = ifelse(status == "carrier", "Carrier", "Non-carrier")) %>%
      rename('Carrier status'='status') %>%
      select(Cohort, `Carrier status`,Gene, Total, Age, BMI, SES, `% Female`, `% Attended University/College`, `Smoking status`, Ethnicity, MDD)
  } else {
  table <- table %>% 
    select(Cohort, Total, Age, BMI, SES, `% Female`, `% Attended University/College`, `Smoking status`, Ethnicity, MDD)
  }
  return(table)
}

demo_both <- rbind(mdd_demo, metabol_demo)
demo_both <- demo_both %>% select(cohort, everything())
View(t(process_demo(demo_both)) %>% as.data.frame())

#### Demographic table separated by carrier status 
demo_tables <- list()
for(i in seq_along(fads_genes$hgnc_symbol)) {
  print(i)
  gene <- fads_genes$hgnc_symbol[i]
  mdd_demo <- read.table(paste0(
    "/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/demographics/",
    gene, "_allcarriers_summarydemo.tsv"), sep= "\t", header =T)
  metabolite_demo <- read.table(paste0(
    "/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/demographics/",
    gene, "_metacarriers_summarydemo.tsv"), sep = "\t", header = T)
  if(gene != "MYRF") {
  mdd_demo <- mdd_demo %>% mutate(cohort = "MDD cohort", Gene = gene) %>% 
    rename('status' = paste0(gene, "_status")) 
  metabolite_demo <- metabolite_demo %>% mutate(cohort ="Metabolite cohort", Gene = gene) %>% 
    rename('status' = paste0(gene, "_status"))
  } else {
    mdd_demo <- mdd_demo %>% mutate(cohort = "MDD cohort", Gene = gene) 
    metabolite_demo <- metabolite_demo %>% mutate(cohort ="Metabolite cohort", Gene = gene) 
  }
 mdd_demo <- process_demo(mdd_demo, TRUE)
 metabolite_demo <- process_demo(metabolite_demo, TRUE)
 demo <- rbind(mdd_demo, metabolite_demo)
 demo <- demo %>% select(Cohort, everything())
demo_tables[[i]] <- demo
}

demo_tables_status <- do.call(rbind, demo_tables)
demo_tables_status <- demo_tables_status %>% select(Gene, Cohort, `Carrier status`, everything())
write.table(demo_tables_status, "/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/demographics/demo_carrier_status.tsv", sep ="\t", row.names =F ,quote =F)
kable(demo_tables_status) %>%
  collapse_rows(columns = c(1,2), valign = "top")




