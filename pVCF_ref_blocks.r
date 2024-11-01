#
.libPaths(c(.libPaths(), "/exports/igmm/eddie/GenScotDepression/users/edavyson/x86_64-pc-linux-gnu-library/4.4"))
library(dplyr)
library(data.table)
library(biomaRt)
library(stringr)

# From https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=2009
# wget on EDDIE and copied to scratch 
# The co-ordinates we are interested in 

# FADS1 gene: 61,799,627-61,829,318
# FADS2 gene: 61,792,980-61,867,354
# FADS3 gene: 61,873,519-61,892,051
# MYRF gene: 61,752,636-61,788,518 
# FEN1 gene: 61,792,911-61,797,238
# TMEM258: 61768501-61792802

# The positions I used in the browser was from the UCSC browser 
# Extract from Ensembl Browser
# Specify the Ensembl dataset you want to use
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand")

# Retrieve gene information
gene_info <- getBM(attributes = attributes, 
                   mart = ensembl)
fads_genes <- gene_info %>% filter(hgnc_symbol %in% c("FADS1", "FADS2", "FADS3", "FEN1", "TMEM258", "MYRF"))
# Transcript information (canonical transcript?)
transcript_attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical")
transcript_info <- getBM(attributes = c(attributes, transcript_attributes), 
                         filters = "hgnc_symbol", 
                         values = fads_genes$hgnc_symbol, 
                         mart = ensembl)
fads_gene_transcript <- transcript_info %>% filter(transcript_is_canonical==1)
gene_positions <- fads_gene_transcript %>% 
mutate(chromosome_name = paste0("chr", chromosome_name)) %>% 
rename(chr=chromosome_name)

write.table(gene_positions, "Resources/FADS_cluster_genes.tsv", sep = "\t", row.names = F, quote = F)

pvcf_ref <- read.table('Resources/dragen_pvcf_coordinates.csv', sep = ',', header = T)
# Calculate the end position of each file (the starting position of the next row -1)
pvcf_ref <- pvcf_ref %>% 
filter(chromosome == "chr11") %>% 
arrange(starting_position) %>% 
mutate(end = lead(starting_position, n = 1, default = NA) -1)

# Filter for the files which start BEFORE the end of the gene 
# OR end AFTER the start of the gene

find_covering_files <- function(gene_start, gene_end, vcf_info) {
    vcf_files <- vcf_info %>% filter((starting_position <= gene_end) & (end >= gene_start)) %>% 
    pull(filename)
    return(vcf_files)
}
# Getting the VCF file(s) for each gene 
gene_coverage <- gene_positions %>% 
rowwise() %>%
mutate(covering_files = list(find_covering_files(start_position, end_position, pvcf_ref))) %>% ungroup()

# Coverting to a character vector
gene_coverage <- gene_coverage %>%
  mutate(covering_files = sapply(covering_files, function(x) paste(x, collapse = ":")))

# Get the numbers of the files 
file_nums <- c(str_split(fads$covering_files, ":")) %>% unlist() %>% unique() %>% sub(".*b(\\d+).*", "\\1", .) %>% as.numeric() %>% sort()

write.table(gene_coverage, "Resources/FADS_cluster_UKB_pVCF.tsv", sep = "\t", row.names = F, quote = F)
# The file ukb24310_c11_b3090_v1.vcf.gz has a starting position of 61799029
# The file ukb24310_c11_b3091_v1.vcf.gz has a starting position of 61819030

# File with information on FADS cluster is b3090_v1 and b3091_v1
