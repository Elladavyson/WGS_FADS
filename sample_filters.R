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

unrelated <- merge(unrelated, missing, by.x = "V1", by.y = "INDV")
print("Unrelated samples of european ancestry with WGS QC'd data for " gene,": ", nrow(unrelated))
# write out the unrelated IDs
print(paste0("Writing out the sample list to unrelated_nomiss_", gene, "_idlist.id"))
readr::write_lines(unrelated$V1, paste0("unrelated_nomiss_", gene, "_idlist.id"))
