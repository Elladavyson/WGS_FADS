## Filtering Quality Control Steps 
# Calculate missingness 

bcftools +fill-tags ${gene}_combined.vcf.gz -Oz -o ${gene}_combined_tag.vcf.gz -- -t F_MISSING

# Filtering variants with 
# - HWE P value < 1e-100
# - Missingness > 0.1

bcftools filter -e 'INFO/HWE < 1e-100 || F_MISSING > 0.1' ${gene}_combined_tag.vcf.gz -o -Oz ${gene}_varqc_combined.vcf.gz 

# Filter to unrelated individuals of european ancestry and to samples < 0.1 missingness 

Rscript sample_filters.R ${gene}

# Outputs a sample list called unrelated_nomiss_${gene}_idlist.id

bcftools view -S unrelated_nomiss_${gene}_idlist.id ${gene}_varqc_combined.vcf.gz -o -Oz ${gene}_varqc_sampqc_comb.vcf.gz

