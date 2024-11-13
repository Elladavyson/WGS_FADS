#!/bin/bash
#########################################################
#$ -N VEP
#$ -l h_rt=48:00:00
#$ -l h_vmem=16G
#$ -cwd
#$ -pe sharedmem 8
#$ -e VEP
#$ -o VEP
#$ -M s2112198@ed.ac.uk
#$ -m baes

. /etc/profile.d/modules.sh
module load igmm/libs/htslib
module load igmm/apps/bcftools
echo 'export PATH=$PATH:/exports/eddie/scratch/s2112198/ensembl-vep' >> ~/.bashrc
source ~/.bashrc
echo "Setting up file paths"
gene=$1
LOF_PATH=/exports/eddie/scratch/s2112198/
ANCESTOR=/exports/eddie/scratch/s2112198/human_ancestor.fa.gz
GERP=/exports/eddie/scratch/s2112198/gerp_conservation_scores.homo_sapiens.GRCh38.bw
CONSVSQL=/exports/eddie/scratch/s2112198/loftee.sql
FASTASEQ=/exports/eddie/scratch/s2112198/GRCh38_full_analysis_set_plus_decoy_hla.fa
CADDSNV=/exports/eddie/scratch/s2112198/whole_genome_SNVs.tsv.gz
CADDINDEL=/exports/eddie/scratch/s2112198/gnomad.genomes.r4.0.indel.tsv.gz # Waiting for this to download from 
REVEL=/exports/eddie/scratch/s2112198/new_tabbed_revel_grch38.tsv.gz
GNOMAD=gnomad.ch.genomesv3.tabbed.tsv.gz

echo "Decompressing the VCF"
bcftools view ${gene}_QC_nosamples.vcf.gz -o ${gene}_QC_nosamples.vcf
echo "Running VEP"
# Get the LoF and CADD 
vep -i ${gene}_QC_nosamples.vcf \
-o annot_vep/${gene}_QC_nosamples_vepannot.vcf \
--vcf \
--canonical \
--fasta $FASTASEQ \
--assembly GRCh38 \
--plugin LoF,loftee_path=$LOF_PATH,human_ancestor_fa=$ANCESTOR,gerp_bigwig=$GERP,conservation_file=$CONSVSQL \
--plugin gnomADc,$GNOMAD \
--plugin REVEL,file=$REVEL \
--custom file=/exports/eddie/scratch/s2112198/gnomad.genomes.v4.1.sites.chr11.vcf.bgz,short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AF_afr%AF_amr%AF_asj%AF_eas%AF_fin%AF_nfe%AF_oth \
--plugin CADD,snv=$CADDSNV,indels=$CADDINDEL \
--fork 4 \
--verbose \
--database \
--force_overwrite \
--assembly GRCh38

echo "Processing VEP output"
# Processing output 
VCF=annot_vep/${gene}_QC_nosamples_vepannot.vcf
 
HEADER=$(printf "chr_pos_ref_alt\nFILTER\n$(bcftools +split-vep ${VCF} -l | cut -f 2)\n" | tr '\n' '\t')

# Splitting the VEP field into different columns 
bcftools +split-vep ${VCF} -d \
-a CSQ \
-A tab \
-f '%CHROM\:%POS\_%REF\_%ALT\t%FILTER\t%CSQ\n'  > ${VCF}_annotated.tsv

echo "$HEADER" | cat - ${VCF}_annotated.tsv > ${VCF}_annotated_hdr.tsv

echo "VEP output in tsv format saved to ${VCF}_annotated_hdr.tsv"

### Downloading files
# Running VEP in EDDIE 
# Git cloning VEP
#git clone https://github.com/Ensembl/ensembl-vep.git
##cd ensembl-vep
# Install on EDDIE
# Loftee Plugin Files here : https://github.com/konradjk/loftee/tree/grch38?tab=readme-ov-file
#echo 'export PATH=$PATH:/exports/eddie/scratch/s2112198/ensembl-vep' >> ~/.bashrc
#source ~/.bashrc

#
# human_ancestor file: https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz
# https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai
# https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi
# GERP file: https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
# FASTA REF FILE : http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
# SQL file (PhyloCSF metrics) : https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz
########### LoFTee files 
#mkdir loftee_hg38/
#cd loftee_hg38/
#wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
#wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz
#wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai
#wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi
#wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz
#gunzip loftee.sql.gz
#cd ../
#tar -czf loftee_hg38.tar.gz loftee_hg38/
#dx upload loftee_hg38.tar.gz

######## GnomAD (for annotation with bcftools annotate)

# You will need to run the following for all chromosomes available, but an example for chromosome 1 is shown below.
# Gnomad is downloaded on EDDIE
# Use the custom flag to get the gnomAD allele frequencies in the same run as CADD, LoFtee and VEP IMPACT 

#/exports/igmm/eddie/BioinformaticsResources/annotation/gnomad/gnomAD_v4.1.0_GRCh38/genome
# GnomAD coverage
#gunzip -c /exports/igmm/eddie/BioinformaticsResources/annotation/gnomad/gnomAD_v3.1.1_GRCh38/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz | sed '1s/.*/#&/' > gnomad.genomesv3.tabbed.tsv
#sed "1s/locus/chr\tpos/; s/:/\t/g" gnomad.genomesv3.tabbed.tsv > gnomad.ch.genomesv3.tabbed.tsv
#bgzip gnomad.ch.genomesv3.tabbed.tsv
#tabix -s 1 -b 2 -e 2 gnomad.ch.genomesv3.tabbed.tsv.gz

##### REVEL
#wget  https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip
#unzip revel-v1.3_all_chromosomes.zip
#cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
#sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
#bgzip new_tabbed_revel.ts
#zcat new_tabbed_revel.tsv.gz | head -n1 > h
#zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
#tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz