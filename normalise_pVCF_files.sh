# Normalise all the pVCF files 

for vcf_file in *.vcf.gz; do
echo "$vcf_file"
base_name=$(basename "$vcf_file" .vcf.gz)
norm_vcf="${base_name}_norm.vcf.gz"
bcftools norm -f GRCh38_full_analysis_set_plus_decoy_hla.fa -m -any -Oz -o "$norm_vcf" "$vcf_file"
bcftools index --tbi "$norm_vcf"
echo "Normalized $vcf_file -> $norm_vcf"
echo "-----------------"
done 


run_normalise="bash normalise_pVCF_files.sh"
 dx run swiss-army-knife \
 -iin="/Bulk/DRAGEN\ WGS/DRAGEN\ population\ level\ WGS\ variants,\ pVCF\ format\ [500k\ release]/chr11/ukb24310_c11_b3091_v1.vcf.gz.tbi" \
 -iin="/Bulk/DRAGEN\ WGS/DRAGEN\ population\ level\ WGS\ variants,\ pVCF\ format\ [500k\ release]/chr11/ukb24310_c11_b3091_v1.vcf.gz" \
 -iin="/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/helper_files/GRCh38_full_analysis_set_plus_decoy_hla.dict" \
 -iin="/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/helper_files/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
 -iin="/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/helper_files/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" \
 -iin="/Code/normalise_pVCF_files.sh" \
 -icmd="${run_normalise}" \
 --tag="normalise_3091" \
 --instance-type "mem1_hdd1_v2_x36" \
 --destination="/Output/" \
 --brief --yes