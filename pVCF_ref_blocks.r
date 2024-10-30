#
library(dplyr)
library(data.table)

# From https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=2009
# wget on EDDIE and copied to scratch 
# The co-ordinates we are interested in 

# FADS1 gene: 61,799,627-61,829,318
# FADS2 gene: 61,792,980-61,867,354
# FADS3 gene: 61,873,519-61,892,051
# MYRF gene: 61,752,636-61,788,518 
# FEN1 gene: 61,792,911-61,797,238

gene_positions <- data.frame(
  Gene = c("FADS1", "FADS2", "FADS3", "MYRF", "FEN1"),
  Start_Position = c(61799627, 61792980, 61873519, 61752636, 61792911),
  End_Position = c(61829318, 61867354, 61892051, 61788518, 61797238)
)

pvcf_ref <- read.table('Resources/dragen_pvcf_coordinates.csv', sep = ',', header = T)
pvcf_chr11 <- pvcf_ref %>% filter(chromosome == 'chr11')
pvcf_chr11 %>% filter(starting_position > 61700000) %>% head()

# The file ukb24310_c11_b3090_v1.vcf.gz has a starting position of 61799029
# The file ukb24310_c11_b3091_v1.vcf.gz has a starting position of 61819030

# File with information on FADS cluster is b3090_v1
