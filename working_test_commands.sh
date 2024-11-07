 Testing it out with FEN1 (only one file)

bcftools view -H ukb24310_c11_b3089_v1.vcf.gz | less
bcftools stats ukb24310_c11_b3089_v1.vcf.gz | less 

# The VCF files have multi-allelic sites and are not left normalised (?)
# Normalisation requires a reference.fa file for left-alignment and normalisation
# https://community.ukbiobank.ac.uk/hc/en-gb/community/posts/17543386269725-Reference-genome-for-DRAGEN-WGS-500k-data
# Which is Bulk / Exome sequences / Exome OQFE CRAM files / helper_files / GRCh38_full_analysis_set_plus_decoy_hla.fa 

# Normalise the pVCF

bcftools view -r chr11:61792911-61797238 ukb24310_c11_b3089_v1.vcf.gz -Ou |
bcftools norm -f GRCh38_full_analysis_set_plus_decoy_hla.fa -m -any -Oz -o ukb24310_c11_b3089_v1_FEN1_norm.vcf.gz 
bcftools index ukb24310_c11_b3089_v1_FEN1_norm.vcf.gz 

# The VCF files have sample information for ALL individuals for every observed rare variant (so ~500,000 rows for each variant)

bcftools view -r chr11:61792911-61797238 ukb24310_c11_b3089_v1.vcf.gz -Ou | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\n' ukb24310_c11_b3089_v1_FEN1_norm.vcf.gz > FEN1_variants.tsv

# Extreracting the genotypes ? 

bcftools query -i 'GT="alt"' -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' ukb24310_c11_b3089_v1_FEN1_norm.vcf.gz > FEN1_var_genotypes.tsv

# Using VEP 

echo -e "CHROM\tPOS\tREF\tALT\t$(bcftools +split-vep -l input.vcf | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > output.tsv
bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\n' -d -A tab ukb24310_c11_b3089_v1_FEN1_norm.vcf.gz >> output.tsv

## Open Cravat
nano dxapp.json # Edit the regionalOptions argument to be "aws:eu-west-2"

# Remove all samples from a VCF
(bcftools query -l ukb24310_c11_b3089_v1_FEN1_norm.vcf.gz) > samples_rm.txt
bcftools view -S ^samples_rm.txt -Oz -o output.vcf.gz ukb24310_c11_b3089_v1_FEN1_norm.vcf.gz
