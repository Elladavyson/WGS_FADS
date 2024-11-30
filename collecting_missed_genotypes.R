library(dplyr)
library(readr)
library(tidyr)

# Download the priority variant list (when filtering to the canonical transcript of the correct gene) e.g FADS2
# dx download /Output/gene_VCF_variants/variants/pri_variants/chrpos/chrposrefalt_gene_canonical/FADS2_priority_annot_chrposrefalt.txt
# Download the priority chr pos file which was used to extract the genotypes (this was not filtered on the FADS2, but rather any canonical transcript of a gene with the more severe consequence)
# Some variants could have been missed
# dx download /Output/gene_VCF_variants/variants/pri_variants/chrpos/FADS2_priority_annot_chrpos.tsv
gene <- "FADS2"
chrpos_file <- read.table(paste0(gene, "_priority_annot_chrpos.tsv"), sep = "\t", header = F)
colnames(chrpos_file) <- c("CHR", "POS")
pri_variants <- readr::read_lines(paste0(gene, "_priority_annot_chrposrefalt.txt"))
pri_variants_df <- data.frame(chr_pos_ref_alt = pri_variants) %>% 
separate(chr_pos_ref_alt, into = c("CHR", "POS", "REF", "ALT"), sep = "_", remove = FALSE)

# Filter to priority variants which are not represented in the chrpos file
# And therefore have not been extracted yet using bcftools
# If there are missing priority variants, save their CHR and POS to extract the genotypes and then append them onto the genotypes file
missed_pri_variants <- pri_variants_df %>% filter(POS %in% chrpos_file$POS==FALSE)
if(nrow(missed_pri_variants) > 0) {
    write.table(missed_pri_variants %>% select(CHR,POS), paste0(gene, "_extra_pri_chrpos.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
}