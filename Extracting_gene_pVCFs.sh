input_file="FADS_genes_remaining_UKB_pVCF.tsv"

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


while IFS=$'\t' read -r hgnc_symbol chr start_position end_position strand ensembl_gene_id ensembl_transcript_id transcript_is_canonical covering_files; do

    # Skip the header line if it's present
    [[ $hgnc_symbol == "hgnc_symbol" ]] && continue
 IFS=":" read -r -a vcf_files <<< "$covering_files"
 region="${chr}:${start_position}-${end_position}"
  # Print the region for the gene
    echo "$hgnc_symbol: $region"
  # Print the VCF files for the gene a
    echo "VCF files for $hgnc_symbol:"
  # For each VCF file 
    for vcf_file in "${vcf_files[@]}"; do
  # Specify file names
    base_name=$(basename "$vcf_file" .vcf.gz)
    norm_file="${base_name}_norm.vcf.gz"
    output_file="${base_name}_norm_${hgnc_symbol}.vcf.gz"
    variant_output_file="${base_name}_${hgnc_symbol}_variants.tsv"
        echo "$vcf_file"
        echo "$norm_file"
  # Filter out the gene region present in the VCF file, normalise and index
        bcftools view -r "$region" "$norm_file" -Oz -o $output_file
        bcftools index --tbi $output_file
        echo "Processed region $region for gene $hgnc_symbol processed from $norm_file and saved to $output_file"
  # Query the variants for Allele Frequency etc 
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\n' $output_file > $variant_output_file
    done
    echo "===================="

 done < "$input_file"

