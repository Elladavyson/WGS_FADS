# Running VEP in EDDIE 
# Git cloning VEP
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
# Get version 110 (fits with the cache files )
cd ensembl-vep
git checkout release/110
perl INSTALL.pl
# Install on EDDIE
# Loftee Plugin Files here : https://github.com/konradjk/loftee/tree/grch38?tab=readme-ov-file
echo 'export PATH=$PATH:/exports/eddie/scratch/s2112198/ensembl-vep' >> ~/.bashrc
source ~/.bashrc

#
# human_ancestor file: https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz
# https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai
# https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi
# GERP file: https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
# FASTA REF FILE : http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
# SQL file (PhyloCSF metrics) : https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz

LOF_PATH=/exports/eddie/scratch/s2112198/
ANCESTOR=/exports/eddie/scratch/s2112198/human_ancestor.fa.gz
GERP=/exports/eddie/scratch/s2112198/gerp_conservation_scores.homo_sapiens.GRCh38.bw
CONSVSQL=/exports/eddie/scratch/s2112198/loftee.sql
FASTASEQ=/exports/eddie/scratch/s2112198/GRCh38_full_analysis_set_plus_decoy_hla.fa
CADDSNV=/exports/eddie/scratch/s2112198/whole_genome_SNVs.tsv.gz
CADDINDEL=/exports/eddie/scratch/s2112198/gnomad.genomes.r4.0.indel.tsv.gz # Waiting for this to download from 
GNOMAD=gnomad.ch.genomesv3.tabbed.tsv.gz
# This is on EDDIE already but you need to run some lines of codes on it (ENSEMBL website)

# Get the LoF and CADD 
vep -i ukb24310_c11_b3089_v1_FEN1_norm_nosamples.vcf \
-o ukb24310_c11_b3089_v1_FEN1_norm_nosamples_vepannot.tsv \
--tab \
--canonical \
--fasta $FASTASEQ \
--plugin LoF,loftee_path=$LOF_PATH,human_ancestor_fa=$ANCESTOR,gerp_bigwig=$GERP,conservation_file=$CONSVSQL \
--plugin gnomADc,$GNOMAD \
--custom file=/exports/eddie/scratch/s2112198/gnomad.genomes.v4.1.sites.chr11.vcf.bgz,short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AF_afr%AF_amr%AF_asj%AF_eas%AF_fin%AF_nfe%AF_oth \
--plugin CADD,snv=$CADDSNV,indels=$CADDINDEL \
--fork 4 \
--verbose \
--database \
--force_overwrite \
--assembly GRCh38


### Downloading files

# LoFTee files 
mkdir loftee_hg38/
cd loftee_hg38/
wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz
wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai
wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi
wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz
gunzip loftee.sql.gz
cd ../
tar -czf loftee_hg38.tar.gz loftee_hg38/
dx upload loftee_hg38.tar.gz

# GnomAD (for annotation with bcftools annotate)

# You will need to run the following for all chromosomes available, but an example for chromosome 1 is shown below.
# Gnomad is downloaded on EDDIE
/exports/igmm/eddie/BioinformaticsResources/annotation/gnomad/gnomAD_v4.1.0_GRCh38/genome
# This is the coverage file 
gunzip -c /exports/igmm/eddie/BioinformaticsResources/annotation/gnomad/gnomAD_v3.1.1_GRCh38/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz | sed '1s/.*/#&/' > gnomad.genomesv3.tabbed.tsv
sed "1s/locus/chr\tpos/; s/:/\t/g" gnomad.genomesv3.tabbed.tsv > gnomad.ch.genomesv3.tabbed.tsv
bgzip gnomad.ch.genomesv3.tabbed.tsv
tabix -s 1 -b 2 -e 2 gnomad.ch.genomesv3.tabbed.tsv.gz

# Use the custom flag to get the gnomAD allele frequencies in the same run as CADD, LoFtee and VEP IMPACT 

