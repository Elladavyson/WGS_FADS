######## RUNNING REGENIE STEP 2 ##########

GENE=$1

### removing the dnanexus paths from the predlist file 

sed 's|/home/dnanexus/out/out/||' metabolite_regenie_step1_pred.list  > basename_metabolite_regenie_step1_pred.list

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
 --phenoFile ukb_unrel_eur_metabol.pheno \
  --covarFile ukb_unrel_eur_covars.covar \
  --phenoCol f.23444.0.0 --phenoCol f.23451.0.0 --phenoCol f.23459.0.0 --phenoCol f.23443.0.0 --phenoCol f.23450.0.0 \
  --covarColList Age,sex_coded,genotype_array,AC,spectrometer,PC{1:10} \
  --catCovarList AC,spectrometer \
  --maxCatLevels 23 \
  --firth --approx \
  --pred basename_metabolite_regenie_step1_pred.list \
  --anno-file annotations_FADS.tsv \
  --set-list setlist_FADS.tsv  \
  --mask-def masks_FADS.txt \
  --aaf-file aaf_FADS.tsv  \
  --aaf-bins 0.01 \
  --set-singletons \
  --extract-setlist ${GENE} \
  --check-burden-files \
  --vc-tests skat,skato,acato \
  --mask-lovo ${GENE},Mask5,0.01 \
  --bsize 100 \
  --out metabolite_${GENE}_step2_LOVO_mask5_0.01

