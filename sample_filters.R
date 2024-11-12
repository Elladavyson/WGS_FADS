if (!requireNamespace("dplyr")) {
  install.packages("dplyr")
}
if (!requireNamespace("optparse")) {
  install.packages("optparse")
}
if (!requireNamespace("readr")) {
  install.packages("readr")
}

library(dplyr)
library(optparse)
library(readr)

parse <- OptionParser()

option_list <- list(
    make_option('--gene', type = "character", help = "Gene name", action = "store")
)
args=commandArgs(trailingOnly = TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args = args)

gene <- opt$gene
# read in the missing file
print(paste0("Reading in the missingness file for the ", gene, " from ", gene, "_missing.miss"))
missing <- read.table(paste0(gene, "_missing.imiss"), header = T)
table(missing$F_MISS)
table(missing$F_MISS < 0.1)
print(paste0("There were ", missing %>% filter(F_MISS > 0.1) %>% nrow(), " individuals with > 0.1 missingness"))

# Read in the list of unrelated individuals of european ancestry

unrelated <- read.table("ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id")

# Get intersection of the samples in the missing and the samples here 

unrelated_missing <- merge(unrelated, missing, by.x = "V1", by.y = "INDV")
print(paste0("Unrelated samples of european ancestry with WGS QC'd data for ", gene,": ", nrow(unrelated)))

unrelated_missing <- unrelated_missing %>% filter(F_MISS <= 0.1)
print(paste0("Unrelated samples of european ancestry and <= 0.1 missingness with WGS QC'd data for ", gene,": ", nrow(unrelated_missing)))
# Also filter to those with a MDD phenotype (PGC3)
mdd <- read.table("MajorDepression.ukb24262.2021-07.txt", header =T)
mdd <- mdd %>% filter(MajDepr != -9)
unrelated_mdd <- merge(unrelated_missing, mdd, by.x = "V1", by.y = "IID")
print(paste0("Unrelated samples of european ancestry with WGS QC'd data for ", gene," and a MDD phenotype (PGC3)", nrow(unrelated_mdd)))

# Write out table of numbers as the sample QC 

sample_QC <- data.table(label= c("Total N", "Unrelated N", "Unrelated N with < 0.1 missingness", "Unrelated N with < 0.1 missingness and MDD phenotype"))
# write out the unrelated IDs
print(paste0("Writing out the sample list to unrelated_nomiss_mdd_", gene, "_idlist.id"))
readr::write_lines(unrelated_mdd$V1, paste0("unrelated_nomiss_mdd_", gene, "_idlist.id"))
