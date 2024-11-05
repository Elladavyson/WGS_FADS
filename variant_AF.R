library(ggalt)
library(stringr)

# dx download *_variants.tsv
# Concatenate all the _variants files together 
system("cat *FADS1_variants.tsv > FADS1_variants.tsv")
system("cat *FADS2_variants.tsv > FADS2_variants.tsv")
system("cat *FADS3_variants.tsv > FADS3_variants.tsv")
system("cat *MYRF_variants.tsv > MYRF_variants.tsv")
system("cat *FEN1_variants.tsv > FEN1_variants.tsv")
system("cat *TMEM258_variants.tsv > TMEM258_variants.tsv")

# Read in the files information 

fads_files <- read.table("FADS_cluster_UKB_pVCF.tsv", sep = "\t", header = T)

# Get the number of the variants in each file (sanity checking results and no truncations)
# FADS1 - 3090 and 3091
fads1_3090 <- read.table("ukb24310_c11_b3090_FADS1_variants.tsv", sep = "\t")
fads1_3091 <- read.table("ukb24310_c11_b3090_FADS1_variants.tsv", sep = "\t")

# FADS2 - 3092 and 3093
fads2_3092 <- read.table("ukb24310_c11_b3092_FADS2_variants.tsv", sep = "\t")
fads2_3093 <- read.table("ukb24310_c11_b3093_FADS2_variants.tsv", sep = "\t")

# FADS3 3093 and 3094
fads3_3093 <- read.table("ukb24310_c11_b3093_FADS3_variants.tsv", sep = "\t")
fads3_3094 <- read.table("ukb24310_c11_b3094_FADS3_variants.tsv", sep = "\t")

# FEN1 3089
fen1_3089 <- read.table("ukb24310_c11_b3089_FEN1_variants.tsv", sep = "\t")

# MYRF 3087, 3088, 3089
myrf_3087 <- read.table("ukb24310_c11_b3087_MYRF_variants.tsv", sep = "\t")
myrf_3088 <- read.table("ukb24310_c11_b3088_MYRF_variants.tsv", sep = "\t")
myrf_3089 <- read.table("ukb24310_c11_b3089_MYRF_variants.tsv", sep = "\t")

# TMEM258 3088 and 3089
tmem_3088 <- read.table("ukb24310_c11_b3088_TMEM258_variants.tsv", sep = "\t")
tmem_3089 <- read.table("ukb24310_c11_b3089_TMEM258_variants.tsv", sep = "\t")
