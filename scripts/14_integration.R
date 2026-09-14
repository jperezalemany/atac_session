library(ggplot2)
library(dplyr)

results_dir <- "results/14_integration"

atac_file <- "results/11_differential/SIM6_vs_CIM7.results.csv"
rna_file <- "resources/SIM6_vs_CIM7.rnaseq.csv"

dir.create(results_dir, showWarnings = F)

message("Reading DARs")
dars <- read.csv(atac_file)
message("Reading DEGs")
degs <- read.csv(rna_file)

df <- merge(dars, degs, by.x = "geneId", by.y = "X", suffixes = c("_ATAC", "_RNA"))

message("DAR-DEG counts:")
# Show counts
df %>%
  filter(dar != "NS" | deg != "NS") %>%
  count(dar, deg)

outfile <- file.path(results_dir, "SIM6_vs_CIM7.atac_rna.csv")
message(paste0("Saving DAR-DEG table to", outfile))
write.csv(df, outfile, quote = F)

# Scatter plot

message("Plotting RNA-ATAC log2 FC scatter plot")
scatter <- df %>%
  filter(dar != "NS" | deg != "NS") %>%
  mutate(dar_deg = paste(dar, deg, sep = "_")) %>%
  ggplot(aes(log2FoldChange_ATAC, log2FoldChange_RNA)) +
  geom_point(aes(color = dar_deg), size = 0.75, alpha = 0.5) +
  labs(x = "SIM6 / CIM7 ATAC", y = "SIM6 / CIM7 RNA") +
  theme_classic()

outfile <- file.path(results_dir, "atac_rna.scatter.png")
message(paste0("Saving scatter plot to ", outfile))
ggsave(outfile, scatter, height = 8, width = 12, units = "cm", dpi = 1200)