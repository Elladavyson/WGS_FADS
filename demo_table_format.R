#### Read in the demographic table 
fads_genes <- read.table("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Intermediate_data_files/FADS_cluster_UKB_pVCF.tsv", sep= "\t", header=T)
mdd_demo <- read.table("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Intermediate_data_files/mdd3_pgc_eur_unrel_sample_demo.tsv", sep = "\t", header = T)
metabol_demo <- read.table("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/Intermediate_data_files/baseline_metabolite_eur_unrel_sample_demo.tsv", sep = "\t", header = T)

mdd_demo <- mdd_demo %>% mutate(cohort = "MDD cohort")
metabol_demo <- metabol_demo %>% mutate(cohort = "Metabolite cohort")

process_demo <- function(table, carrier_status, meta = FALSE) {
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
                              'Total' = 'total',
                              'Ethnicity'= 'ethnicity'
                            )
  if(meta == TRUE) {
    print("Changing metabolite names")
    table <- table %>%
      rename("Degree of Unsaturation"="f.23443.0.0_stat",
             "Docosahexaenoic Acid" = "f.23450.0.0_stat",
             "Omega-3 Fatty Acids" ="f.23444.0.0_stat",
             "Omega-3 Fatty Acids to Total Fatty Acids %"="f.23451.0.0_stat",
             "Omega-6 Fatty Acids to Omega-3 Fatty Acids Ratio" = "f.23459.0.0_stat"
      )
  }
  if(carrier_status == "Yes") {
    table <- table %>% 
      mutate(status = ifelse(status == "carrier", "Carrier", "Non-carrier")) %>%
      rename('Carrier status'='status',
             'Cohort' = 'cohort') %>%
      select(Cohort, `Carrier status`,Gene, Total, Age, BMI, SES, `% Female`, `% Attended University/College`, `Smoking status`, Ethnicity, MDD)
  } else if(carrier_status == "No") {
  table <- table %>% 
    rename('Cohort' = 'cohort')
    select(Cohort, Total, Age, BMI, SES, `% Female`, `% Attended University/College`, `Smoking status`, Ethnicity, MDD)
  } else if(carrier_status == "lovo"){
    if(meta == TRUE) {
      table <- table %>%
        select(ends_with("carrier"), Total, Age, BMI, SES, `% Female`, `% Attended University/College`, `Smoking status`, Ethnicity, MDD,
               `Degree of Unsaturation`, `Docosahexaenoic Acid`, `Omega-3 Fatty Acids`,
               `Omega-3 Fatty Acids to Total Fatty Acids %`, `Omega-6 Fatty Acids to Omega-3 Fatty Acids Ratio`)
    } else {
    table <- table %>% 
      select(ends_with("carrier"), Total, Age, BMI, SES, `% Female`, `% Attended University/College`, `Smoking status`, Ethnicity, MDD)
}
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
 mdd_demo <- process_demo(mdd_demo, "Yes")
 metabolite_demo <- process_demo(metabolite_demo, "Yes")
 demo <- rbind(mdd_demo, metabolite_demo)
 demo <- demo %>% select(Cohort, everything())
demo_tables[[i]] <- demo
}

demo_tables_status <- do.call(rbind, demo_tables)
demo_tables_status <- demo_tables_status %>% select(Gene, Cohort, `Carrier status`, everything())
write.table(demo_tables_status, "/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/demographics/demo_carrier_status.tsv", sep ="\t", row.names =F ,quote =F)
kable(demo_tables_status) %>%
  collapse_rows(columns = c(1,2), valign = "top")

### LOVO carriers 

lovo_variants <- c("11_61816814_G_C", "11_61810815_C_A")

mdd_lovo_demo <- read.table("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/demographics/FADS1_chr11_61816814_G_C_demo.tsv", sep = "\t", header=T)
mdd_lovo_demo <- mdd_lovo_demo %>% 
  rename('chr11_61816814_G_C carrier'='lovo_carrier')
View(t(process_demo(mdd_lovo_demo, "lovo")) %>% as.data.frame())

metabolite_lovo_demo <- read.table("/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/demographics/FADS1_chr11_61810815_C_A_meta_demo.tsv", sep = "\t", header=T)
metabolite_lovo_demo <- metabolite_lovo_demo %>% 
  rename('chr11_61810815_C_A carrier'='lovo_carrier')
metabolite_lovo_demo <- t(process_demo(metabolite_lovo_demo, "lovo", meta = TRUE)) %>% 
  as.data.frame() %>% 
  mutate(name = rownames(.)) %>% 
  select(name, everything())
colnames(metabolite_lovo_demo) <- metabolite_lovo_demo[1,]
metabolite_lovo_demo <- metabolite_lovo_demo[-1, ]
  
write.table(metabolite_lovo_demo, 
            "/Users/ellad/UniversityEdinburgh/PhD/Year_3/WGS_proj/demographics/chr11_61810815_C_A_meta_proc_demo.tsv", sep = "\t", quote = F, row.names = F)

