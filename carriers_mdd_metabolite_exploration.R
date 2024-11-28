library(data.table)
library(dplyr)

# carriers info 

FEN1_carriers <- read.table("FEN1_carrier_info.tsv", sep = "\t", header = T)
mdd <- read.table("MajorDepression.ukb24262.2021-07.txt", header = T)

# Merge with MDD

FEN1_carriers_mdd <- merge(FEN1_carriers, mdd, by.x = "SAMPLE", by.y = "IID")