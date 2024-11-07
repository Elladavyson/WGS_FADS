
### SAK submission commands 
## Normalising 

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

### Extracting for each gene 

run_extract="bash Extracting_gene_pVCFs.sh {genename}"
 dx run swiss-army-knife \
 -iin="/Output/ukb24310_c11_b3087_v1_norm.vcf.gz.tbi" \
 -iin="/Output/ukb24310_c11_b3087_v1_norm.vcf.gz" \
  -iin="/Output/ukb24310_c11_b3088_v1_norm.vcf.gz.tbi" \
 -iin="/Output/ukb24310_c11_b3088_v1_norm.vcf.gz" \
   -iin="/Output/ukb24310_c11_b3089_v1_norm.vcf.gz.tbi" \
 -iin="/Output/ukb24310_c11_b3089_v1_norm.vcf.gz" \
 -iin="{genename}_UKB_pVCF.tsv" \
 -iin="/Code/Extracting_gene_pVCFs.sh"  \
 -icmd="${run_extract}" \
 --tag="extract_{genename}_3087_3088_3089" \
 --instance-type "mem1_hdd1_v2_x36" \
 --destination="/Output/" \
 --brief --yes


run_extract="bash Extracting_gene_pVCFs.sh TMEM258"
 dx run swiss-army-knife \
  -iin="/Output/ukb24310_c11_b3088_v1_norm.vcf.gz.tbi" \
 -iin="/Output/ukb24310_c11_b3088_v1_norm.vcf.gz" \
   -iin="/Output/ukb24310_c11_b3089_v1_norm.vcf.gz.tbi" \
 -iin="/Output/ukb24310_c11_b3089_v1_norm.vcf.gz" \
 -iin="TMEM258_UKB_pVCF.tsv" \
 -iin="/Code/Extracting_gene_pVCFs.sh"  \
 -icmd="${run_extract}" \
 --tag="extract_TMEM258_3087_3088_3089" \
 --instance-type "mem1_hdd1_v2_x36" \
 --destination="/Output/" \
 --brief --yes

run_extract="bash Extracting_gene_pVCFs.sh FEN1"
 dx run swiss-army-knife \
   -iin="/Output/ukb24310_c11_b3089_v1_norm.vcf.gz.tbi" \
 -iin="/Output/ukb24310_c11_b3089_v1_norm.vcf.gz" \
 -iin="FEN1_UKB_pVCF.tsv" \
 -iin="/Code/Extracting_gene_pVCFs.sh"  \
 -icmd="${run_extract}" \
 --tag="extract_FEN1_3087_3088_3089" \
 --instance-type "mem1_hdd1_v2_x36" \
 --destination="/Output/" \
 --brief --yes

 run_extract="bash Extracting_gene_pVCFs.sh FADS3"
 dx run swiss-army-knife \
   -iin="/Output/ukb24310_c11_b3093_v1_norm.vcf.gz.tbi" \
 -iin="/Output/ukb24310_c11_b3093_v1_norm.vcf.gz" \
  -iin="/Output/ukb24310_c11_b3094_v1_norm.vcf.gz.tbi" \
 -iin="/Output/ukb24310_c11_b3094_v1_norm.vcf.gz" \
 -iin="FADS3_UKB_pVCF.tsv" \
 -iin="/Code/Extracting_gene_pVCFs.sh"  \
 -icmd="${run_extract}" \
 --tag="extract_FADS3_3087_3088_3089" \
 --instance-type "mem1_hdd1_v2_x36" \
 --destination="/Output/" \
 --brief --yes
