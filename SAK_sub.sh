
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
 --instance-type "mem1_hdd1_v2_x72" \
 --destination="/Output/" \
 --brief --yes

### Extracting for each gene 

run_extract="bash Extracting_gene_pVCFs.sh MYRF"
 dx run swiss-army-knife \
 -iin="/Output/ukb24310_c11_b3087_v1_norm.vcf.gz.tbi" \
 -iin="/Output/ukb24310_c11_b3087_v1_norm.vcf.gz" \
  -iin="/Output/ukb24310_c11_b3088_v1_norm.vcf.gz.tbi" \
 -iin="/Output/ukb24310_c11_b3088_v1_norm.vcf.gz" \
   -iin="/Output/ukb24310_c11_b3089_v1_norm.vcf.gz.tbi" \
 -iin="/Output/ukb24310_c11_b3089_v1_norm.vcf.gz" \
 -iin="MYRF_UKB_pVCF.tsv" \
 -iin="/Code/Extracting_gene_pVCFs.sh"  \
 -icmd="${run_extract}" \
 --tag="extract_MYRF_3087_3088_3089" \
 --instance-type "mem1_hdd1_v2_x36" \
 --destination="/Output/" \
 --brief --yes


run_extract="bash Extracting_gene_pVCFs.sh FADS1"
 dx run swiss-army-knife \
 -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3089_v1_norm.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3089_v1_norm.vcf.gz" \
  -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3090_v1_norm.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3090_v1_norm.vcf.gz" \
   -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3091_v1_norm.vcf.gz.tbi" \
   -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3091_v1_norm.vcf.gz" \
 -iin="/Input/FADS1_UKB_pVCF.tsv" \
 -iin="/Code/Extracting_gene_pVCFs.sh"  \
 -icmd="${run_extract}" \
 --tag="extract_FADS1_3089_3090_3091" \
 --instance-type "mem1_hdd1_v2_x36" \
 --destination="/Output/" \
 --brief --yes

 run_extract="bash Extracting_gene_pVCFs.sh FADS2"
 dx run swiss-army-knife \
 -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3089_v1_norm.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3089_v1_norm.vcf.gz" \
  -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3090_v1_norm.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3090_v1_norm.vcf.gz" \
   -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3091_v1_norm.vcf.gz.tbi" \
   -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3091_v1_norm.vcf.gz" \
      -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3092_v1_norm.vcf.gz.tbi" \
   -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3092_v1_norm.vcf.gz" \
      -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3093_v1_norm.vcf.gz.tbi" \
   -iin="/Output/gene_VCF_variants/normalised_vcf/ukb24310_c11_b3093_v1_norm.vcf.gz" \
 -iin="/Input/FADS2_UKB_pVCF.tsv" \
 -iin="/Code/Extracting_gene_pVCFs.sh"  \
 -icmd="${run_extract}" \
 --tag="extract_FADS2_3089_3090_3091_3092_3093" \
 --instance-type "mem1_hdd1_v2_x72" \
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
 --instance-type "mem1_hdd1_v2_x72" \
 --destination="/Output/" \
 --brief --yes

## Concatenating the FADS1 and FADS2 gene VCFs 
 run_concat="bcftools concat -o FADS2_combined.vcf.gz -O z ukb24310_c11_b3089_v1_norm_FADS2.vcf.gz ukb24310_c11_b3090_v1_norm_FADS2.vcf.gz ukb24310_c11_b3091_v1_norm_FADS2.vcf.gz ukb24310_c11_b3092_v1_norm_FADS2.vcf.gz ukb24310_c11_b3093_v1_norm_FADS2.vcf.gz"
 dx run swiss-army-knife \
   -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3089_v1_norm_FADS2.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3089_v1_norm_FADS2.vcf.gz" \
  -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3090_v1_norm_FADS2.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3090_v1_norm_FADS2.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3091_v1_norm_FADS2.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3091_v1_norm_FADS2.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3092_v1_norm_FADS2.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3092_v1_norm_FADS2.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3093_v1_norm_FADS2.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3093_v1_norm_FADS2.vcf.gz" \
 -icmd="${run_concat}" \
 --tag="concat_FADS2_3089_3090_3091_3092_3093" \
 --instance-type "mem1_hdd1_v2_x72" \
 --destination="/Output/" \
 --brief --yes

 run_concat="bcftools concat -o FADS1_combined.vcf.gz -O z  ukb24310_c11_b3090_v1_norm_FADS2.vcf.gz ukb24310_c11_b3091_v1_norm_FADS2.vcf.gz"
 dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3090_v1_norm_FADS2.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3090_v1_norm_FADS2.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3091_v1_norm_FADS2.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3091_v1_norm_FADS2.vcf.gz" \
 -icmd="${run_concat}" \
 --tag="concat_FADS1_3090_3091" \
 --instance-type "mem1_hdd1_v2_x72" \
 --destination="/Output/" \
 --brief --yes

## Indexing and extracting format fields 

run_index="bcftools index --tbi FADS2_combined.vcf.gz"
 dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FADS2_combined.vcf.gz" \
 -icmd="${run_index}" \
 --tag="index_FADS2" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/gene_VCF_variants/gene_vcfs/" \
 --brief --yes

 run_index="bcftools index --tbi FADS1_combined.vcf.gz"
 dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FADS1_combined.vcf.gz" \
 -icmd="${run_index}" \
 --tag="index_FADS1" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/gene_VCF_variants/gene_vcfs/" \
 --brief --yes

# Getting the sample missingness information from vcftools 
 run_sample_missingness="vcftools --gzvcf FEN1_combined.vcf.gz --missing-indv --out FEN1_missing"
  dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FEN1_combined.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/FEN1_combined.vcf.gz.tbi" \
 -icmd="${run_sample_missingness}" \
 --tag="missingness_FEN1" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes

  run_sample_missingness="vcftools --gzvcf FADS1_combined.vcf.gz --missing-indv --out FADS1_missing"
  dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FADS1_combined.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/FADS1_combined.vcf.gz.tbi" \
 -icmd="${run_sample_missingness}" \
 --tag="missingness_FADS1" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes
   run_sample_missingness="vcftools --gzvcf FADS3_combined.vcf.gz --missing-indv --out FADS3_missing"
  dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FADS3_combined.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/FADS3_combined.vcf.gz.tbi" \
 -icmd="${run_sample_missingness}" \
 --tag="missingness_FADS3" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes

    run_sample_missingness="vcftools --gzvcf MYRF_combined.vcf.gz --missing-indv --out MYRF_missing"
  dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/MYRF_combined.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/MYRF_combined.vcf.gz.tbi" \
 -icmd="${run_sample_missingness}" \
 --tag="missingness_MYRF" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes

     run_sample_missingness="vcftools --gzvcf TMEM258_combined.vcf.gz --missing-indv --out TMEM258_missing"
  dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/TMEM258_combined.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/TMEM258_combined.vcf.gz.tbi" \
 -icmd="${run_sample_missingness}" \
 --tag="missingness_TMEM258" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes

 ## Running variant level and sample level QC

 run_qc="bash qc_filters.sh FEN1"
   dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FEN1_combined.vcf.gz" \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FEN1_combined.vcf.gz.tbi" \
  -iin="/Output/gene_VCF_variants/QualControl/FEN1_missing.imiss" \
  -iin="/Input/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id" \
  -iin="/Input/MajorDepression.ukb24262.2021-07.txt" \
  -iin="/Code/qc_filters.sh" \
  -iin="/Code/sample_filters.R" \
 -icmd="${run_qc}" \
 --tag="FEN1_QC" \
 --instance-type "mem1_hdd1_v2_x36" \
 --destination="/Output/" \
 --brief --yes