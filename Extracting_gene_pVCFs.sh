input_file="test_cluster_UKB_pVCF.tsv"

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
  # Get the base name for file tracking (i.e the field ID in UKB)
    base_name=$(basename "$vcf_file" .vcf.gz)
  # Define the output file names (one VCF and one TSV)
    output_file="${basename}_norm_${hgnc_symbol}.vcf.gz"
    variant_output_file="${basename}_${hgnc_symbol}_variants.tsv"
        echo "$vcf_file"
  # Filter out the gene region present in the VCF file, normalise and index
        bcftools view -r "$region" "$vcf_file" -Ou |
        bcftools norm -f GRCh38_full_analysis_set_plus_decoy_hla.fa -m -any -Oz -o $output_file
        bcftools index --tbi $output_file
        echo "Normalised, left-aligned and biallelic sites for $hgnc_symbol processed from $vcf_file and saved to $output_file"
  # Query the variants for Allele Frequency etc 
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\n' $output_file > $variant_output_file
        dx upload --path "$DX_PROJECT_CONTEXT_ID:" $variant_output_file
    done
    echo "===================="

 done < "$input_file"

