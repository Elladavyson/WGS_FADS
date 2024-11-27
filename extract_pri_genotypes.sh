## Extracting prioritised variants genotypes from QC VCF file

gene="${1}"

## Index the file 

bcftools index --tbi ${gene}_varqc_sampqc_HWE_combined.vcf.gz

## Extract the genotypes of the prioritised variants 
bcftools view -R ${gene}_priority_annot_chrpos.tsv ${gene}_varqc_sampqc_HWE_combined.vcf.gz | bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' > ${gene}_priority_genotypes.tsv