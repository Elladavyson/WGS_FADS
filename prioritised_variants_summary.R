library(tidyr)
library(dplyr)

# Read in the VEP output for the genes 

genes <- c("FADS1", "FADS2", "FADS3", "FEN1", "MYRF", "TMEM258")
for (i in genes) {
    print(i)
vepoutput <- read.table(paste0("AC_overzero/", i, "_vepoutput_proc/", i, "_priority_annot.tsv"),
sep = '\t', header = T)
vepoutput$GENE = i
assign(paste0(i, "_vep_pri"), vepoutput)
}

# Create a function to calculate counts in each AF range
count_af_ranges <- function(df) {
  range1 <- nrow(df[df$AF < 0.000001, ])        
    range2 <- nrow(df[df$AF >= 0.000001 & df$AF < 0.00001, ])   
  range3 <- nrow(df[df$AF >= 0.00001 & df$AF < 0.0001, ])                         # AF < 0.0001
  range4 <- nrow(df[df$AF >= 0.0001 & df$AF < 0.001, ])               # 0.0001 ≤ AF < 0.001
  range5 <- nrow(df[df$AF >= 0.001 & df$AF < 0.01, ])                  # 0.001 ≤ AF < 0.1
  return(c(range1, range2, range3, range4, range5))
}

count_ac_ranges <- function(df) {
    range1 <- nrow(df[df$AC == 1, ])
    range2 <- nrow(df[df$AC > 1 & df$AC < 5, ])
    range3 <- nrow(df[df$AC >= 5 & df$AC < 10, ])
    range4 <- nrow(df[df$AC >=10 & df$AC < 50, ])
    range5 <- nrow(df[df$AC >= 50 & df$AC < 100, ])
    range6 <- nrow(df[df$AC >= 100 & df$AC < 1000, ])
    range7 <- nrow(df[df$AC >=1000, ])
    return(c(range1, range2, range3, range4, range5, range6, range7))
}

# Apply this function to each gene's data frame     
af_ranges <- sapply(list(FADS1_vep_pri, FADS2_vep_pri, FADS3_vep_pri, FEN1_vep_pri, 
                      MYRF_vep_pri, TMEM258_vep_pri), count_af_ranges)
ac_ranges <- sapply(list(FADS1_vep_pri, FADS2_vep_pri, FADS3_vep_pri, FEN1_vep_pri, 
                      MYRF_vep_pri, TMEM258_vep_pri), count_ac_ranges)
# Create a data frame with results
af_counts <- data.frame(
    Gene = genes,
    Range1_AF_less_0.000001 = af_ranges[1, ],
    Range2_AF_0.000001_to_0.00001 = af_ranges[2, ],
    Range3_AF_0.00001_to_0.0001 = af_ranges[3, ],
    Range4_AF_0.0001_to_0.001 = af_ranges[4, ],
    Range5_AF_0.001_to_0.01 = af_ranges[5, ]
)
ac_counts <- data.frame(
    Gene = genes,
    Range1_AC_1 = ac_ranges[1, ],
    Range2_AC_1_5= ac_ranges[2, ],
    Range3_AC_5_10 = ac_ranges[3, ],
    Range4_AC_10_50 = ac_ranges[4, ],
    Range5_AC_50_100 = ac_ranges[5, ],
    Range6_AC_100_1000 = ac_ranges[6, ],
    Range7_AC_over_1000 = ac_ranges[7, ]
)
# View the final table
print(af_counts)
print(ac_counts)

n_variants <- data.frame(Gene = genes, All = c(nrow(FADS1_vep_pri), nrow(FADS2_vep_pri), nrow(FADS3_vep_pri), 
nrow(FEN1_vep_pri), nrow(MYRF_vep_pri), nrow(TMEM258_vep_pri)))

n_variants <- merge(n_variants, af_counts, by = "Gene")
n_variants <- merge(n_variants, ac_counts, by = "Gene")
long_data <- pivot_longer(
  n_variants,
  cols = c(starts_with("Range"), "All"),  # Select all columns starting with "Range"
  names_to = "AF_AC_Range",        # New column for the range labels
  values_to = "Variant_Count"   # New column for the counts
)

# View the resulting data
print(long_data)

## Plot 

plot_af <- ggplot(long_data %>% filter(grepl("AF", AF_AC_Range)), aes(x = AF_AC_Range, y = Variant_Count, fill = AF_AC_Range)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Variant_Count), 
            vjust = -0.5,           # Adjust text position above the bar
            size = 3.5) +          # Text size
  facet_wrap(~Gene) +
  theme_minimal() +
  labs(
    title = "Variant Counts by Allele Frequency Range for Each Gene",
    x = "Allele Frequency Range",
    y = "Variant Count"
  ) +
  theme(
    axis.text.x = element_blank(), # Remove x-axis tick labels
    axis.ticks.x = element_blank(), # Remove x-axis ticks (optional)
    axis.text.x.bottom = element_blank() # Removes from bottom x-axis only (optional)
  )


plot_ac <- ggplot(long_data %>% filter(grepl("AC", AF_AC_Range)), aes(x = AF_AC_Range, y = Variant_Count, fill = AF_AC_Range)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Variant_Count), 
            vjust = -0.5,           # Adjust text position above the bar
            size = 3.5) +          # Text size
  facet_wrap(~Gene) +
  theme_minimal() +
  labs(
    title = "Variant Counts by Allele Count Range for Each Gene",
    x = "Allele Count Range",
    y = "Variant Count"
  ) +
  theme(
    axis.text.x = element_blank(), # Remove x-axis tick labels
    axis.ticks.x = element_blank(), # Remove x-axis ticks (optional)
    axis.text.x.bottom = element_blank() # Removes from bottom x-axis only (optional)
  )

png("AF_gene_priority_variants.png", 
    width = 3000, height = 2000, res = 300, type = "cairo")
plot_af
dev.off()

png("AC_gene_priority_variants.png",
width = 3000, height = 2000, res = 300, type = "cairo")
plot_ac
dev.off()