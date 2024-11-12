## Check the high missingness reported in FADS1 by vcftools 

# Randomly sample 100 individuals with a reported missing ness > 0.1 by VCF tools 
FADS1_missing <- read.table("FADS1_missing.imiss", header =T)
sample_fads_missing <- sample_n(FADS1_missing, 100)
readr::write_lines(sample_fads_missing$INDV, "sample_missingFADS1.id")

# Extract the genotypes for these individuals from the combined gene VCF (no QC steps applied)
# I.e the same file provided to VCFtools to report on missingness

# query_genotypes="bcftools view -S sample_missingFADS1.id FADS1_combined.vcf.gz | bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' > FADS1_nsample_oQC_missing_genotypes.tsv"
 #  dx run swiss-army-knife \
 # -iin="/Output/gene_VCF_variants/gene_vcfs/FADS1_combined.vcf.gz" \
 # -iin="/Output/gene_VCF_variants/gene_vcfs/FADS1_combined.vcf.gz.tbi" \
 # -iin="sample_missingFADS1.id" \
## -icmd="${query_genotypes}" \
 #--tag="FADS1_geno_missingness_check" \
 #--instance-type "mem1_hdd1_v2_x36" \
 #--destination="/Output/" \
 #--brief --yes

 # Read the output back into R

FADS1_samples <- read.table("FADS1_nsample_oQC_missing_genotypes.tsv", sep = "\t")
colnames(FADS1_samples) <- c("CHR", "POS", "REF", "ALT", "SAMPLE", "GT")
FADS1_samples$chr_pos_ref_alt <- paste0(FADS1_samples$CHR, "_", FADS1_samples$POS, "_", FADS1_samples$REF, "_", FADS1_samples$ALT)
# FADS1 samples missingness 
FADS1_samples_missing <- FADS1_samples %>% 
group_by(SAMPLE) %>% 
summarise(n_missing = sum(GT =="." | GT == "./.")) %>% 
mutate(prop_missing = n_missing/13548)

# Cross reference these values back with those reported from VCFtools
FADS1_missing_both <- merge(FADS1_missing, FADS1_samples_missing, by.x = 'INDV', by.y = "SAMPLE")

table(FADS1_missing_both$F_MISS== round(FADS1_missing_both$prop_missing,6))
FADS1_sample_missingness <- ggplot(FADS1_missing_both, aes(x = F_MISS, y= prop_missing)) + 
geom_point() + 
labs(x = "VCFTools Missingness", y = "Manual Missingness" , title= "Test of missingness in FADS1 in 100 samples") + 
theme_minimal()

ggsave(filename = "FADS1_missingness_sample.png", FADS1_sample_missingness, width = 6, height = 6, device = "png", dpi = 300)

#dx upload FADS1_missingness_sample.png
