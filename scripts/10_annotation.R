#!/usr/bin/env Rscript

library(ChIPseeker)
library(txdbmaker)
library(rtracklayer)

library(dplyr)
library(ggplot2)
library(tidyr)

gtf <- "resources/Araport11_GTF_genes_transposons.20241001.subset.gtf"
summits <- "results/09_counts/all_peaks.summits.bed"
tss <- c(-500, 500)
results_dir <- "results/10_annotation"

dir.create(results_dir, showWarnings = F)

txdb <- makeTxDbFromGFF(gtf)

ann <- annotatePeak(
  summits, tssRegion = tss, TxDb = txdb, level = "gene", overlap = "all"
) %>%
  as.data.frame() %>%
  mutate(feature = sub(" \\([A-Za-z0-9 ,./]+\\)", "", annotation)) %>%
  mutate(distanceToTES = if_else(geneStrand == 1, end - geneEnd, geneStart - start)) %>%
  mutate(min_region = if_else(abs(distanceToTSS) < abs(distanceToTES), "TSS", "TES")) %>%
  mutate(feature_simp = case_when(
    min_region == "TSS" & abs(distanceToTSS) <= 500 ~ "TSS",
    min_region == "TES" & abs(distanceToTES) <= 500 ~ "TES",
    distanceToTSS > 0 & distanceToTES < 0 ~ "Gene Body",
    distanceToTSS < -500 | distanceToTES > 500 ~ "Intergenic"
  ))


# Make pie chart
message("Plotting pie chart of genomic features")
pie_chart <- ann %>%
  group_by(feature) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = "", y = n, fill = feature)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white")) +
  theme(legend.title = element_blank())

# Save pie chart
message(paste0("Saving pie chart to ", file.path(results_dir, "peaks_pie_chart.png")))
ggsave(
  file.path(results_dir, "peaks_pie_chart.png"),
  pie_chart,
  height = 8, width = 12, units = "cm", dpi = 1200
)

# Plot simplified pie chart
pie_chart2 <- ann %>%
  group_by(feature_simp) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = "", y = n, fill = feature_simp)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white")) +
  theme(legend.title = element_blank())

# Save pie chart
message(paste0("Saving pie chart to ", file.path(results_dir, "peaks_pie_chart_simple.png")))
ggsave(
  file.path(results_dir, "peaks_pie_chart_simple.png"),
  pie_chart2,
  height = 8, width = 12, units = "cm", dpi = 1200
)

# Save table
message(paste0("Saving annotated peaks to ", file.path(results_dir, "all_peaks.annotation.txt")))
write.table(
  ann, file.path(results_dir, "all_peaks.annotation.txt"),
  sep = "\t", row.names = F, quote = F
)
