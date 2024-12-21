curl -O -L https://github.com/rgcgithub/regenie/releases/download/v3.2.6/regenie_v3.2.6.gz_x86_64_Linux_mkl.zip; unzip regenie_v3.2.6.gz_x86_64_Linux_mkl.zip -d /usr/local/bin/; chmod a+x /usr/local/bin/regenie_v3.2.6.gz_x86_64_Linux_mkl; rm regenie*.zip
regenie_v3.2.6.gz_x86_64_Linux_mkl --step 2 \
 --bed ${GENE}_varqc_sampqc_HWE_combined \
  --phenoFile ukb_unrel_eur_pgc3_mdd.pheno \
  --covarFile ukb_unrel_eur_covars.covar \
  --phenoCol MajDepr \
  --covarColList Age,sex_coded,genotype_array,AC,spectrometer,PC{1:10} \
  --catCovarList AC,spectrometer \
  --maxCatLevels 23 \
  --bt \
  --firth --approx \
  --pred basename_mdd_regenie_step1_BT_pred.list \
  --compute-corr \
  --output-corr-text \
  --bsize 20 \
  --out mdd_${GENE}_step2_BT_LD_cor