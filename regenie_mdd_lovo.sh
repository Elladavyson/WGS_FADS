######## RUNNING REGENIE STEP 2 ##########

GENE=$1

### removing the dnanexus paths from the predlist file 

sed 's|/home/dnanexus/out/out/||' mdd_regenie_step1_pred.list > basename_mdd_regenie_step1_pred.list

## Using the updated bim files
# Manually in R 
#
#bim <- read.table("{gene}_varqc_sampqc_HWE_combined.bim", header =F)
#colnames(bim) <- c("CHR", "ID", "DIST", "POS", "A1", "A2")
#bim$ID = paste0(bim$CHR, ":", bim$POS, ":", bim$A2, ":", bim$A1)
#write.table(bim, "{gene}_varqc_sampqc_HWE_combined.bim", row.names =F, quote = F, col.names = F)

# Downloading latest version of REGENIE
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
  --anno-file annotations_FADS.tsv \
  --set-list setlist_FADS.tsv  \
  --mask-def masks_FADS.txt \
  --aaf-file aaf_FADS.tsv  \
  --aaf-bins 0.01 \
  --set-singletons \
  --extract-setlist ${GENE} \
  --check-burden-files \
  --vc-tests skat,skato,acato \
  --mask-lovo ${GENE},Mask4,0.01 \
  --bsize 100 \
  --out mdd_${GENE}_step2_LOVO