## Filtering Quality Control Steps 
# Calculate missingness 

echo "Getting the F_missing tag per variant"
bcftools +fill-tags ${gene}_combined.vcf.gz -Oz -o ${gene}_combined_tag.vcf.gz -- -t F_MISSING

# Filtering variants with 
# - HWE P value < 1e-100
# - Missingness > 0.1
echo "Filtering on the HWE < 1E-100 and F_MISISNG > 0.1"
bcftools filter -e 'INFO/HWE < 1e-100 || F_MISSING > 0.1' ${gene}_combined_tag.vcf.gz -o -Oz ${gene}_varqc_combined.vcf.gz 

# Filter to unrelated individuals of european ancestry and to samples < 0.1 missingness 
echo "Filtering to unrelated individuals with sample missingness < 0.1 and of european ancestry"
Rscript sample_filters.R ${gene}

# Outputs a sample list called unrelated_nomiss_${gene}_idlist.id

bcftools view -S unrelated_nomiss_${gene}_idlist.id -o ${gene}_varqc_combined.vcf.gz -Oz ${gene}_varqc_sampqc_comb.vcf.gz

echo "Filtered and QC variants in ${gene}_varqc_sampqc_comb.vcf.gz"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\n' ${gene}_varqc_sampqc_comb.vcf.gz > ${gene}_variants_qc.vcf.gz

echo "Filtered and QC variants extracted and saved in tabular format at ${gene}_variants_qc.vcf.gz"