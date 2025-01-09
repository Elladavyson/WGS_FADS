
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
  -iin="/Output/gene_VCF_variants/gene_vcfs/raw_geneVCFs/FADS2_combined.vcf.gz" \
  -iin="/Output/gene_VCF_variants/gene_vcfs/raw_geneVCFs/FADS2_combined.vcf.gz.tbi" \
  -iin="/Output/gene_VCF_variants/QualControl/sample_missingness/FADS2_missing.imiss" \
  -iin="/Input/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id" \
  -iin="/Input/MajorDepression.ukb24262.2021-07.txt" \
  -iin="/Code/qc_filters.sh" \
  -iin="/Code/sample_filters.R" \
 -icmd="${run_qc}" \
 --tag="FADS2_QC" \
 --instance-type "mem1_hdd1_v2_x36" \
 --destination="/Output/" \
 --brief --yes

 run_query="bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\n' FADS1_QC_nosamples.vcf.gz > FADS1_QC_variants.tsv"
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

  run_HWEfilter="bcftools filter -e 'INFO/HWE < 1e-100 || INFO/HWE == 0 || F_MISSING > 0.1' MYRF_varqc_sampqc_comb.vcf.gz -Oz -o MYRF_varqc_sampqc_HWE_combined.vcf.gz"
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

 extract_qc_metrics="bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/F_MISSING\t%INFO/HWE\n' FADS2_combined_tag.vcf.gz > FADS2_all_qcmetrics.tsv"
dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/QC/FADS2_combined_tag.vcf.gz" \
 -icmd="${extract_qc_metrics}" \
 --tag="QC metrics" \
 --instance-type "mem1_ssd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes

  strip_samples="bcftools view -S ^unrelated_nomiss_mdd_FADS2_idlist.id FADS2_varqc_sampqc_comb.vcf.gz -Oz -o FADS2_QC_nosamples.vcf.gz"
dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/QC/FADS2_varqc_sampqc_comb.vcf.gz" \
  -iin="/Output/gene_VCF_variants/QualControl/sample_list/unrelated_nomiss_mdd_FADS2_idlist.id" \
 -icmd="${strip_samples}" \
 --tag="FADS2 strip samples" \
 --instance-type "mem1_hdd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes

 check_zeroAC="bcftools query -R FEN1_QC_var_zeroAC.tsv -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n' FEN1_combined.vcf.gz > FEN1_AC_all_qczero.tsv"
     dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/gene_vcfs/raw_geneVCFs/FEN1_combined.vcf.gz" \
  -iin="/Output/gene_VCF_variants/gene_vcfs/raw_geneVCFs/FEN1_combined.vcf.gz.tbi" \
  -iin="FEN1_QC_var_zeroAC.tsv" \
 -icmd="${check_zeroAC}" \
 --tag="FEN1_AC_check_for_variants_withACzero_inQC" \
 --instance-type "mem1_ssd1_v2_x8" \
 --destination="/Output/" \
 --brief --yes

query_pri_genotypes="bash extract_pri_genotypes.sh FADS3"
  dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/variants/pri_variants/FADS3_priority_annot_chrpos.tsv" \
  -iin="/Output/gene_VCF_variants/gene_vcfs/QC/FADS3_varqc_sampqc_HWE_combined.vcf.gz" \
  -iin="/Code/extract_pri_genotypes.sh" \
  -icmd="${query_pri_genotypes}" \
 --tag="Extract_FADS3_priority_genotypes" \
 --instance-type "mem1_ssd1_v2_x4" \
 --destination="/Output/" \
 --brief --yes

 query_extra_pri_genotypes="bash extra_genotypes_extract.sh TMEM258"
  dx run swiss-army-knife \
  -iin="/Output/gene_VCF_variants/variants/pri_variants/chrpos/TMEM258_extra_pri_chrpos.tsv" \
  -iin="/Output/gene_VCF_variants/gene_vcfs/QC/TMEM258_varqc_sampqc_HWE_combined.vcf.gz" \
  -iin="/Code/extra_genotypes_extract.sh" \
  -icmd="${query_extra_pri_genotypes}" \
 --tag="Extract_extra_TMEM258_priority_genotypes" \
 --instance-type "mem1_ssd1_v2_x8" \
 --destination="/Output/genotypes/" \
 --brief --yes


 summarise_genotypes="Rscript priority_carriers_counts.R --gene FADS3 --extra Yes"
   dx run swiss-army-knife \
  -iin="/Code/priority_carriers_counts.R" \
  -iin="/Input/FADS_cluster_UKB_pVCF.tsv" \
  -iin="/Output/gene_VCF_variants/variants/pri_variants/chrpos/chrposrefalt_gene_canonical/FADS3_priority_annot_chrposrefalt.txt" \
  -iin="/Output/genotypes/FADS3_priority_genotypes.tsv" \
  -iin="/Output/genotypes/FADS3_extra_priority_genotypes.tsv" \
  -icmd="${summarise_genotypes}" \
 --tag="Summarise_FADS3_genotypes" \
 --instance-type "mem1_ssd1_v2_x16" \
 --destination="/Output/genotypes/genotype_summary/" \
 --brief --yes

 merge_genotype_files="cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]* .;\
ls *.bed | sed -e 's/.bed//g'> files_to_merge.txt;\
plink --merge-list files_to_merge.txt --make-bed \
--autosome-xy --out ukb22418_c1_22_v2_merged;\
rm files_to_merge.txt;"

dx run swiss-army-knife \
  -icmd="${merge_genotype_files}" \
 --tag="merge genotype calls" \
 --instance-type "mem1_ssd1_v2_x16" \
 --destination="/Output/regenie/input/" \
 --brief --yes

 run_plink_qc="plink2 --bfile ukb22418_c1_22_v2_merged\
 --keep ukb_unrel_eur_metabol.pheno --autosome\
 --maf 0.01 --mac 20 --geno 0.1\
 --hwe 1e-15 --mind 0.1\
 --write-snplist --write-samples\
 --no-id-header --out \
 WGS_array_snps_qc_pass"

 dx run swiss-army-knife \
  -icmd="${run_plink_qc}" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.bed" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.bim" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.fam" \
  -iin="/Output/regenie/ukb_unrel_eur_metabol.pheno" \
 --tag="qc merged genotypes" \
 --instance-type "mem1_ssd1_v2_x16" \
 --destination="/Output/regenie/input/" \
 --brief --yes

 regenie_step1_meta="regenie --step 1\
 --lowmem --out metabolite_regenie_step1 --bed ukb22418_c1_22_v2_merged\
 --phenoFile ukb_unrel_eur_metabol.pheno --covarFile ukb_unrel_eur_covars.covar\
 --extract WGS_array_snps_qc_pass.snplist\
 --phenoCol f.23444.0.0 --phenoCol f.23451.0.0 --phenoCol f.23459.0.0 --phenoCol f.23443.0.0 --phenoCol f.23450.0.0\
 --covarColList Age,sex_coded,genotype_array,AC,spectrometer,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10\
 --catCovarList AC,spectrometer\
 --maxCatLevels 23\
 --bsize 1000 --loocv --gz --threads 16"

  dx run swiss-army-knife \
  -icmd="${regenie_step1_meta}" \
  -iin="/Output/regenie/input/WGS_array_snps_qc_pass.snplist" \
  -iin="/Output/regenie/ukb_unrel_eur_metabol.pheno" \
  -iin="/Output/regenie/ukb_unrel_eur_covars.covar" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.bed" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.bim" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.fam" \
 --tag="regenie step 1 metabolite" \
 --instance-type "mem1_ssd1_v2_x16" \
 --destination="/Output/regenie/input/" \
 --brief --yes

convert_vcf_bed="plink2 --vcf FEN1_varqc_sampqc_HWE_combined.vcf.gz --make-bed -out FEN1_varqc_sampqc_HWE_combined"
dx run swiss-army-knife \
  -icmd="${convert_vcf_bed}" \
  -iin="/Output/gene_VCF_variants/gene_vcfs/QC/FEN1_varqc_sampqc_HWE_combined.vcf.gz" \
  -iin="/Output/gene_VCF_variants/gene_vcfs/QC/FEN1_varqc_sampqc_HWE_combined.vcf.gz.tbi" \
 --tag="convert FEN1 VCF" \
 --instance-type "mem1_ssd1_v2_x16" \
 --destination="/Output/regenie/input/WGS_BEDS/plink2" \
 --brief --yes

####### Running REGENIE step 2 using SAK (REGENIE VERSION 3.1.1)

dx run swiss-army-knife \
  -icmd="bash regenie_step2.sh FADS2" \
  -iin="/Output/regenie/input/WGS_BEDS/FADS2_varqc_sampqc_HWE_combined.bed" \
  -iin="/Output/regenie/input/WGS_BEDS/FADS2_varqc_sampqc_HWE_combined.bim" \
  -iin="/Output/regenie/input/WGS_BEDS/FADS2_varqc_sampqc_HWE_combined.fam" \
  -iin="/Output/regenie/ukb_unrel_eur_metabol.pheno" \
  -iin="/Output/regenie/ukb_unrel_eur_covars.covar" \
  -iin="/Output/regenie/input/annotations_FADS.tsv" \
  -iin="/Output/regenie/input/masks_FADS.txt" \
  -iin="/Output/regenie/input/aaf_FADS.tsv" \
  -iin="/Output/regenie/input/setlist_FADS.tsv" \
  -iin="/Output/regenie/input/metabolite_regenie_step1_1.loco.gz" \
  -iin="/Output/regenie/input/metabolite_regenie_step1_2.loco.gz" \
  -iin="/Output/regenie/input/metabolite_regenie_step1_3.loco.gz" \
  -iin="/Output/regenie/input/metabolite_regenie_step1_4.loco.gz" \
  -iin="/Output/regenie/input/metabolite_regenie_step1_5.loco.gz" \
  -iin="/Output/regenie/input/metabolite_regenie_step1_pred.list" \
  -iin="/Code/regenie_step2.sh" \
  --tag="regenie step 2 FADS2" \
  --instance-type "mem1_ssd1_v2_x16" \
  --destination="/Output/regenie/step2_res/version3.2.6/metabolite/minMAC/" \
  --brief --yes

######## REGENIE QC of MDD 

 run_plink_qc_mdd="plink2 --bfile ukb22418_c1_22_v2_merged\
 --keep ukb_unrel_eur_pgc3_mdd.pheno --autosome\
 --maf 0.01 --mac 20 --geno 0.1\
 --hwe 1e-15 --mind 0.1\
 --write-snplist --write-samples\
 --no-id-header --out \
 WGS_array_snps_qc_pass_pgc3_mdd"

 dx run swiss-army-knife \
  -icmd="${run_plink_qc_mdd}" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.bed" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.bim" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.fam" \
  -iin="/Output/regenie/ukb_unrel_eur_pgc3_mdd.pheno" \
 --tag="qc merged genotypes MDD phenotype" \
 --instance-type "mem1_ssd1_v2_x16" \
 --destination="/Output/regenie/input/" \
 --brief --yes

## REGENIE STEP 1 MDD

 regenie_step1_mdd="regenie --step 1\
 --lowmem --out mdd_regenie_step1_BT_nospectrometer --bed ukb22418_c1_22_v2_merged\
 --phenoFile ukb_unrel_eur_pgc3_mdd.pheno --covarFile ukb_unrel_eur_covars.covar\
 --extract WGS_array_snps_qc_pass_pgc3_mdd.snplist\
 --phenoCol MajDepr \
 --covarColList Age,sex_coded,genotype_array,AC,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --catCovarList AC \
 --maxCatLevels 23\
 --bt \
 --firth --approx \
 --bsize 1000 --loocv --gz --threads 16"

  dx run swiss-army-knife \
  -icmd="${regenie_step1_mdd}" \
  -iin="/Output/regenie/input/WGS_array_snps_qc_pass_pgc3_mdd.snplist" \
  -iin="/Output/regenie/ukb_unrel_eur_pgc3_mdd.pheno" \
  -iin="/Output/regenie/ukb_unrel_eur_covars.covar" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.bed" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.bim" \
  -iin="/Output/regenie/input/GRCh37_merged/ukb22418_c1_22_v2_merged.fam" \
 --tag="regenie step 1 mdd" \
 --instance-type "mem1_ssd1_v2_x16" \
 --destination="/Output/regenie/input/" \
 --brief --yes

##### REGENIE STEP 2 MDD


dx run swiss-army-knife \
  -icmd="bash regenie_mdd_step2.sh TMEM258" \
  -iin="/Output/regenie/input/WGS_BEDS/TMEM258_varqc_sampqc_HWE_combined.bed" \
  -iin="/Output/regenie/input/WGS_BEDS/TMEM258_varqc_sampqc_HWE_combined.bim" \
  -iin="/Output/regenie/input/WGS_BEDS/TMEM258_varqc_sampqc_HWE_combined.fam" \
  -iin="/Output/regenie/ukb_unrel_eur_pgc3_mdd.pheno" \
  -iin="/Output/regenie/ukb_unrel_eur_covars.covar" \
  -iin="/Output/regenie/input/annotations_FADS.tsv" \
  -iin="/Output/regenie/input/masks_FADS.txt" \
  -iin="/Output/regenie/input/aaf_FADS.tsv" \
  -iin="/Output/regenie/input/setlist_FADS.tsv" \
  -iin="/Output/regenie/input/mdd_regenie_step1_BT_nospectrometer_1.loco.gz" \
  -iin="/Output/regenie/input/mdd_regenie_step1_BT_nospectrometer_pred.list" \
  -iin="/Code/regenie_mdd_step2.sh" \
  --tag="regenie step 2 MDD TMEM258" \
  --instance-type "mem1_ssd1_v2_x16" \
  --destination="/Output/regenie/step2_res/version3.2.6/MDD/bt_pred_nospectrometer/" \
  --brief --yes

# LOVO run FOR mask 4 0.01

dx run swiss-army-knife \
  -icmd="bash regenie_mdd_lovo.sh FADS1" \
  -iin="/Output/regenie/input/WGS_BEDS/FADS1_varqc_sampqc_HWE_combined.bed" \
  -iin="/Output/regenie/input/WGS_BEDS/FADS1_varqc_sampqc_HWE_combined.bim" \
  -iin="/Output/regenie/input/WGS_BEDS/FADS1_varqc_sampqc_HWE_combined.fam" \
  -iin="/Output/regenie/ukb_unrel_eur_pgc3_mdd.pheno" \
  -iin="/Output/regenie/ukb_unrel_eur_covars.covar" \
  -iin="/Output/regenie/input/annotations_FADS.tsv" \
  -iin="/Output/regenie/input/masks_FADS.txt" \
  -iin="/Output/regenie/input/aaf_FADS.tsv" \
  -iin="/Output/regenie/input/setlist_FADS.tsv" \
  -iin="/Output/regenie/input/mdd_regenie_step1_BT_nospectrometer_1.loco.gz" \
  -iin="/Output/regenie/input/mdd_regenie_step1_BT_nospectrometer_pred.list" \
  -iin="/Code/regenie_mdd_lovo.sh" \
  --tag="regenie step 2 MDD LOVO FADS1" \
  --instance-type "mem1_ssd1_v2_x16" \
  --destination="/Output/regenie/step2_res/version3.2.6/MDD/lovo/" \
  --brief --yes

 ## LD heatmaps

dx run swiss-army-knife \
  -icmd="Rscript LD_heatmaps.R" \
  -iin="/Output/regenie/regenie_LD/mdd_FADS1_step2_BT_LD_cor.corr" \
  -iin="/Output/regenie/regenie_LD/mdd_FADS1_step2_BT_LD_cor.corr.snplist" \
  -iin="/Output/regenie/regenie_LD/mdd_FADS2_step2_BT_LD_cor.corr" \
  -iin="/Output/regenie/regenie_LD/mdd_FADS2_step2_BT_LD_cor.corr.snplist" \
  -iin="/Output/regenie/regenie_LD/mdd_FADS3_step2_BT_LD_cor.corr" \
  -iin="/Output/regenie/regenie_LD/mdd_FADS3_step2_BT_LD_cor.corr.snplist" \
  -iin="/Output/regenie/regenie_LD/mdd_MYRF_step2_BT_LD_cor.corr" \
  -iin="/Output/regenie/regenie_LD/mdd_MYRF_step2_BT_LD_cor.corr.snplist" \
  -iin="/Output/regenie/regenie_LD/mdd_FEN1_step2_BT_LD_cor.corr" \
  -iin="/Output/regenie/regenie_LD/mdd_FEN1_step2_BT_LD_cor.corr.snplist" \
  -iin="/Output/regenie/regenie_LD/mdd_TMEM258_step2_BT_LD_cor.corr" \
  -iin="/Output/regenie/regenie_LD/mdd_TMEM258_step2_BT_LD_cor.corr.snplist" \
  -iin="/Input/FADS_cluster_UKB_pVCF.tsv" \
  -iin="/Output/gene_VCF_variants/variants/pri_variants/chrpos/chrposrefalt_gene_canonical/FADS1_priority_annot_chrposrefalt.txt" \
  -iin="/Output/gene_VCF_variants/variants/pri_variants/chrpos/chrposrefalt_gene_canonical/FADS2_priority_annot_chrposrefalt.txt" \
  -iin="/Output/gene_VCF_variants/variants/pri_variants/chrpos/chrposrefalt_gene_canonical/FADS3_priority_annot_chrposrefalt.txt" \
  -iin="/Output/gene_VCF_variants/variants/pri_variants/chrpos/chrposrefalt_gene_canonical/FEN1_priority_annot_chrposrefalt.txt" \
  -iin="/Output/gene_VCF_variants/variants/pri_variants/chrpos/chrposrefalt_gene_canonical/MYRF_priority_annot_chrposrefalt.txt" \
  -iin="/Output/gene_VCF_variants/variants/pri_variants/chrpos/chrposrefalt_gene_canonical/TMEM258_priority_annot_chrposrefalt.txt" \
  -iin="/Code/LD_heatmaps.R" \
  --tag="LD heatmaps" \
  --instance-type "mem1_ssd1_v2_x16" \
  --destination="/Output/regenie/regenie_LD/" \
  --brief --yes
