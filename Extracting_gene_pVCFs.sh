input_file="${1}_UKB_pVCF.tsv"
echo $input_file 
while IFS=$'\t' read -r hgnc_symbol chr start_position end_position strand ensembl_gene_id ensembl_transcript_id transcript_is_canonical covering_files; do

    # Skip the header line if it's present
    if [[ "$hgnc_symbol" == "hgnc_symbol" ]]; then 
    continue
    fi
 IFS=":" read -r -a vcf_files <<< "$covering_files"
 region="${chr}:${start_position}-${end_position}"
  # Print the region for the gene
    echo "$hgnc_symbol: $region"
  # Print the VCF files for the gene a
    echo "VCF files for $hgnc_symbol:"

    # Temporary list to hold the output VCF files (indexed for gene region)

    gene_vcf_list=()
  # For each VCF file 
    for vcf_file in "${vcf_files[@]}"; do
  # Specify file names
    base_name=$(basename "$vcf_file" .vcf.gz)
    norm_file="${base_name}_norm.vcf.gz"
    output_file="${base_name}_norm_${hgnc_symbol}.vcf.gz"
        echo "Processing $vcf_file for region $region"
        echo "Normalised file: $norm_file"
  # Filter out the gene region present in the VCF file, normalise and index
        bcftools view -r "$region" "$norm_file" -Oz -o $output_file
        bcftools index --tbi $output_file
        echo "Processed region $region for gene $hgnc_symbol processed from $norm_file and saved to $output_file"
        gene_vcf_list+=("$output_file")
    done 
  # Query the variants for Allele Frequency etc 
  # Concatenate all gene-specific indexed VCFs
  concatenated_vcf="${hgnc_symbol}_combined.vcf.gz"
  bcftools concat -Oz -o "$concatenated_vcf" "${gene_vcf_list[@]}"
  bcftools index --tbi "$concatenated_vcf"

  echo "Concatenated VCF for $hgnc_symbol saved as $concatenated_vcf"
  variant_output_file="${hgnc_symbol}_variants.tsv"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\n' $concatenated_vcf > $variant_output_file
  echo "Queried variants for gene $hgnc_symbol saved to $variant_output_file"
  echo "===================="

 done < "$input_file"

