
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

 run_concat="bcftools concat -o FADS1_combined.vcf.gz -O z  ukb24310_c11_b3090_v1_norm_FADS1.vcf.gz ukb24310_c11_b3091_v1_norm_FADS1.vcf.gz"
 dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3090_v1_norm_FADS1.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3090_v1_norm_FADS1.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3091_v1_norm_FADS1.vcf.gz.tbi" \
 -iin="/Output/gene_VCF_variants/gene_vcfs/ukb24310_c11_b3091_v1_norm_FADS1.vcf.gz" \
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

     run_sample_missingness="vcftools --gzvcf FADS1_combined.vcf.gz --missing-indv --out FADS1_missing"
  dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FADS1_combined.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/FADS1_combined.vcf.gz.tbi" \
 -icmd="${run_sample_missingness}" \
 --tag="missingness_FADS1" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/gene_VCF_variants/QualControl/" \
 --brief --yes

 ## Running variant level and sample level QC

 run_qc="bash qc_filters.sh FADS2"
   dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FADS2_combined.vcf.gz" \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FADS2_combined.vcf.gz.tbi" \
  -iin="/Output/gene_VCF_variants/QualControl/FADS2_missing.imiss" \
  -iin="/Input/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id" \
  -iin="/Input/MajorDepression.ukb24262.2021-07.txt" \
  -iin="/Code/qc_filters.sh" \
  -iin="/Code/sample_filters.R" \
 -icmd="${run_qc}" \
 --tag="FADS2_QC" \
 --instance-type "mem1_hdd1_v2_x36" \
 --destination="/Output/" \
 --brief --yes

 run_query="bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\n' FADS2_combined.vcf.gz > FADS2_variants.tsv"
 dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FADS2_combined.vcf.gz" \
   -iin="/Output/gene_VCF_variants/gene_vcfs/FADS2_combined.vcf.gz.tbi" \
 -icmd="${run_query}" \
 --tag="query_FADS2" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/gene_VCF_variants/gene_vcfs/" \
 --brief --yes

# The QC script I began first did not use the correct id list file, so rederived them and then manually doing the sample QC step again

  run_sampqc="bcftools view -S unrelated_nomiss_mdd_TMEM258_idlist.id TMEM258_varqc_combined.vcf.gz -Oz -o TMEM258_varqc_sampqc_comb.vcf.gz"
 dx run swiss-army-knife \
  -iin="/Output/TMEM258_varqc_combined.vcf.gz" \
  -iin="/Output/unrelated_nomiss_mdd_TMEM258_idlist.id" \
 -icmd="${run_sampqc}" \
 --tag="sample_filter_TMEM258" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes

# The QC script doesn't seem to be filtering on the HWE flag properly, so doing this for each gene too 

  run_HWEfilter="bcftools filter -e 'INFO/HWE < 1e-100 || INFO/HWE == 0 || F_MISSING > 0.1' FEN1_varqc_sampqc_comb.vcf.gz -Oz -o FEN1_varqc_sampqc_HWE_combined.vcf.gz"
 dx run swiss-army-knife \
  -iin="/Output/FEN1_varqc_sampqc_comb.vcf.gz" \
 -icmd="${run_HWEfilter}" \
 --tag="HWE_filter_FEN1" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes

# There is high missingness in the FADS1 gene from the VCFTools input- manually checking this missingness in the genotypes! 

 query_genotypes="bcftools view -S sample_missingFADS1.id FADS1_combined.vcf.gz | bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' > FADS1_nsample_oQC_missing_genotypes.tsv"
    dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FADS1_combined.vcf.gz" \
  -iin="/Output/gene_VCF_variants/gene_vcfs/FADS1_combined.vcf.gz.tbi" \
  -iin="sample_missingFADS1.id" \
 -icmd="${query_genotypes}" \
 --tag="FADS1_geno_missingness_check" \
 --instance-type "mem1_hdd1_v2_x36" \
 --destination="/Output/" \
 --brief --yes

 extract_qc_metrics="bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/F_MISSING\t%INFO/HWE\n" MYRF_combined_tag.vcf.gz > MYRF_all_qcmetrics.tsv"
dx run swiss-army-knife \
  -iin="/Output/MYRF_combined_tag.vcf.gz" \
 -icmd="${extract_qc_metrics}" \
 --tag="QC metrics" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes

  strip_samples="bcftools view -S ^unrelated_nomiss_mdd_FEN1_idlist.id FEN1_varqc_sampqc_HWE_combined.vcf.gz -Oz -o FEN1_QC_nosamples.vcf.gz"
dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/QC/FEN1_varqc_sampqc_HWE_combined.vcf.gz" \
  -iin="/Output/gene_VCF_variants/QualControl/sample_list/unrelated_nomiss_mdd_FEN1_idlist.id" \
 -icmd="${strip_samples}" \
 --tag="FEN1 strip samples" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes