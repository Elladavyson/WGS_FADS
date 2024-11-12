# dx download the variants 
# dx download /Input/FADS_cluster_UKB_pVCF.tsv
# dx download /Output/gene_VCF_variants/QualControl/sample_missingness/*
# dx download /Output/gene_VCF_variants/QualControl/variant_QC/*
library(ggplot2)
library(dplyr)
# Read in the FADS coordinates used to index the VCF files 
fads <- read.table("FADS_cluster_UKB_pVCF.tsv", sep = "\t", header = T)
for (i in c(1:nrow(fads))) {
    gene <- fads$hgnc_symbol[i]
    print(gene)
    print(paste0("gene_co-ordinates: ", fads$start_position[i], "-", fads$end_position[i]))
print(paste0("Reading in the missingness file for the ", gene, " from ", gene, "_missing.miss"))
missing <- read.table(paste0(gene, "_missing.imiss"), header = T)
table(missing$F_MISS)
table(missing$F_MISS < 0.1)
print(paste0("There were ", missing %>% filter(F_MISS > 0.1) %>% nrow(), " individuals with > 0.1 missingness"))

# Read in the list of unrelated individuals of european ancestry

unrelated <- read.table("ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.id")

# Get intersection of the samples in the missing and the samples here 

unrelated_missing <- merge(unrelated, missing, by.x = "V1", by.y = "INDV")
missingness_hist <- ggplot(unrelated_missing, aes(x = F_MISS)) + 
geom_histogram() + geom_vline(xintercept  = 0.1, color = "red", linetype = "dashed")+ 
theme_minimal() + labs(x = "Proportion of variants missing", y = "Count", title = paste0(gene, ": Unrelated EUR, N = ", nrow(unrelated_missing)))
assign(paste0(gene, "_missing_hist"), missingness_hist)
}
ggsave(filename = "missingness_hists.png", ggarrange(MYRF_missing_hist, FEN1_missing_hist, FADS1_missing_hist, FADS2_missing_hist, FADS3_missing_hist, TMEM258_missing_hist), width = 10, height = 8, device = "png", dpi = 300)

tmem_missing <- read.table("TMEM258_sampleQC_num.tsv", sep = '\t', header =T )
fads1_missing <- read.table("FADS1_sampleQC_num.tsv", sep = '\t', header =T )
fads2_missing <- read.table("FADS2_sampleQC_num.tsv", sep = '\t', header = T)
myrf_missing <- read.table("MYRF_sampleQC_num.tsv", sep = '\t', header =T )
fen1_missing <- read.table("FEN1_sampleQC_num.tsv", sep = '\t', header =T )
fads3_missing <- read.table("FADS3_sampleQC_num.tsv", sep = '\t', header =T )
all_sample_qc <- rbind(fads1_missing, fads2_missing, fads3_missing, tmem_missing, fen1_missing, myrf_missing)
write.table(all_sample_qc, "all_fads_sample_qc.tsv", sep = "\t", row.names = F, quote = F)

sample_qc_long <-  all_sample_qc %>%
    pivot_longer(cols = -GENE, 
                 names_to = "Category", 
                 values_to = "Count")

sample_qc_bar <- ggplot(sample_qc_long, aes(x = GENE, y = Count, fill = Category)) + 
geom_bar(stat="identity", position=position_dodge(width = 0.8)) +  
geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.7, size = 2) + 
theme_minimal() + 
theme(legend.position = "top")
ggsave(filename = "sample_qc_summarybar.png", sample_qc_bar, width = 8, height = 6, device = "png", dpi = 300)

# dx upload sample_qc_summarybar.png
# dx upload missingness_hists.png
# dx upload all_sample_qc.tsv
# dx upload all_fads_sample_qc.tsv
fads <- fads %>% filter(hgnc_symbol != "FADS2")

for (i in c(1:nrow(fads))) {
    gene <- fads$hgnc_symbol[i]
    print(gene)
    print(paste0("gene_co-ordinates: ", fads$start_position[i], "-", fads$end_position[i]))
    variant_qc <- read.table(paste0(gene, "_all_qcmetrics.tsv"), sep = "\t", header =T)
    colnames(variant_qc) <- c("CHR", "POS", "REF", "ALT", "F_MISSING", "HWE")
    variant_qc_summary <- variant_qc %>%
    summarize(
        QC_passed =sum(F_MISSING <= 0.1 & HWE > 1E-100),
        F_missing_gt_0.1 = sum(F_MISSING > 0.1),
    HWE_lt_1e_100= sum(HWE < 1e-100),
    both_HWE_F_MISSING = sum(HWE < 1e-100 & F_MISSING > 0.1)
    ) %>%
     mutate(GENE = gene) %>% 
     select(GENE, everything())
 assign(paste0(gene, "_qc_summary"), variant_qc_summary)

}

variant_qc_summary <- rbind(FADS1_qc_summary,FADS3_qc_summary, TMEM258_qc_summary, MYRF_qc_summary, FEN1_qc_summary)

variant_qc_long <-  variant_qc_summary %>%
    pivot_longer(cols = -GENE, 
                 names_to = "Category", 
                 values_to = "Count")

variant_qc_bar <- ggplot(variant_qc_long, aes(x = GENE, y = Count, fill = Category)) + 
geom_bar(stat="identity", position = "dodge") +  
geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.7) + 
theme_minimal() + 
theme(legend.position = "top") +
scale_fill_discrete(
labels = c("QC_passed"="QC passed",
"F_missing_gt_0.1" = "F MISSING > 0.1",
"HWE_lt_1e_100" = "HWE < 1e-100",
"both_HWE_F_MISSING" = "F_MISSING > 0.1 & HWE <1E-100"))

ggsave(filename = "FADS_variant_qc_bar.png", variant_qc_bar , width = 10, height = 6, device = "png", dpi = 300)
# dx upload FADS_variant_qc_bar.png