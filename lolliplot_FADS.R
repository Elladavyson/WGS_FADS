################# Create a lolliplot of prioritised variants ###################
library(BiocManager)
library(data.table)
library(dplyr)

if (!requireNamespace("trackViewer", quietly = TRUE)) {
  BiocManager::install("trackViewer")
}
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  BiocManager::install("rtracklayer")
}
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  BiocManager::install("GenomicRanges")
}
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
library(trackViewer)
library(optparse)

#### Read in the data 
Ensembl_file <- "ExonsSpreadsheet-Homo_sapiens_Transcript_Exons_ENST00000305885.csv"
ensembl <- read_csv(Ensembl_file) %>% 
as.data.frame() %>% 
select(-Sequence) %>% 
filter(!grepl('sequence', `Exon / Intron`)) %>%
mutate(ExonIntron = ifelse(grepl('Intron',`Exon / Intron`), "Intron", "Exon"))%>% 
mutate(color_for_lolliplot = ifelse(ExonIntron=="Intron","#51C6E6","#FF8833"))

# Carrier genotype info
carriers_genotypes <- read.table(paste0(gene, "_carriers_genotypes.tsv"), sep = '\t', header = T)
# Variant carrier summary
variant_carrier_summary <- read.table(paste0(gene, "_variant_count_summary.tsv"), sep = "\t", header = T)
strand_str = "Forward"

######### Making genomic ranges object for gene 
features <- GRanges("chr11", IRanges(start=c(ensembl$Start),
end=c(ensembl$End),
names = c(ensembl$ExonIntron),
))
features$fill <-  c(ensembl$color_for_lolliplot)

######### Make genomic ranges object for the prioritised variants 

# Establish the compound heterozygous variants 
poten_comp_het_variants <- carriers_genotypes %>% 
filter(n_het_variants > 1 & n_althom_variants == 0 ) %>% 
pull(het_variants) %>% strsplit(., ':') %>%
unlist() %>% 
unique()

variant_carrier_summary <- variant_carrier_summary %>% 
mutate(zygosity = case_when(chr_pos_ref_alt %in% poten_comp_het_variants ~ "Inconc_comp_het",
althet_carriers > 0 & althom_carriers > 0 ~ "Alt_hom_and_alt_het",
althet_carriers > 0 & althom_carriers ==0 ~ "Alt-het",
althet_carriers == 0 & althom_carriers > 0 ~ "Alt-hom",
TRUE ~ "No valid zygosity"))

variant <- variant_carrier_summary %>% mutate(colour_for_lolliplot = 
case_when(zygosity == "Alt-het" ~ "#90EE90",
zygosity == "Alt_hom_and_alt_het" ~ "#006400", 
zygosity == "Alt-hom" ~ "#ADD8E6",
zygosity == "Inconc_comp_het" ~  "lightpink",
TRUE ~ "No valid zygosity"))

######### Ranges object for the variants 
SNP <- as.numeric(variant$POS)
sample.gr <- GRanges(unique(variant$CHR), 
IRanges(SNP, width=1, names=paste0(SNP,"_",variant$REF, "_", variant$ALT)), 
strand=ifelse(strand_str == "Forward", "+", "-"))
sample.gr$color <- variant$colour_for_lolliplot

print("Variant Granges object")
sample.gr

######### Color names for the legend 
color_dict <- c("Alt-het" = "#90EE90",
"Alt_hom_and_alt_het" = "#006400",
"Alt-hom" = "#ADD8E6", 
"Inconc_comp_het" = "lightpink")

get_color_name <- function(color_codes) {
    color_names <- sapply(color_codes, function(color_code) {
    pos <- which(color_dict == color_code)
    if (length(pos) > 0) {
        return(names(color_dict)[pos])
    } else {
        return(NA)
    }
})
return(color_names)
}
legend <- sample.gr$color %>% unique()# legend fill color
names(legend) <- c(get_color_name(legend))# legend labels
legend <- list(labels=c("Alt-het" , "Inconc_comp_het"), 
col = c("#90EE90","lightpink"), 
fill = c("#90EE90","lightpink"))

######### Plotting it 

scores <- variant_carrier_summary %>% 
arrange(factor(chr_pos_ref_alt %in% paste0("chr11_", names(sample.gr)))) %>% 
pull(althet_carriers)
# Get the number of carriers 
names(sample.gr)
sample.gr$score <- scores +100
sample.gr$node.label <- as.character(sample.gr$score - 100)
sample.gr.rot <- sample.gr
sample.gr.rot$label.parameter.label <- ""
sample.gr.rot$shape <- "circle"

## The legend and the y axis label are not showing ? 

print(paste0("Saving the lolliplot to ", outdir, "lolliplot_variants_all.png"))
png(paste0(outdir, "lolliplot_variants_all.png"),
    width = 3600, height = 1200, res = 300, type = "cairo")
lolliplot(sample.gr.rot, features, legend=legend, cex = 0.001, yaxis = FALSE, ylab="Number of Carriers in UKB Sample")
grid.text("Prioritised variants in FEN1", x=.5, y=.9, just="top",
          gp=gpar(cex=1, fontface="bold"))
dev.off()

