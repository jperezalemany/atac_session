
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

message("\n==== 07_tss_enrichment.R ====\n")

results_dir <- "results/07_tss_enrichment"
dir.create(results_dir, showWarnings = F)

read_deeptools_profile <- function(filename) {
  data <- read.table(filename, skip = 2)
  total_bins <- ncol(data) - 2
  colnames(data) <- c("sample", "genes", c(1:total_bins))
  
  data %>%
    pivot_longer(as.character(c(1:total_bins)), names_to = "bin", values_to = "signal") %>%
    mutate(bin = as.integer(bin))
}

estimate_background <- function(tss) {
  total_bins <- max(tss$bin)
  binsize <- 2000 / total_bins
  flank <- 100 / binsize
  
  tss %>%
    group_by(sample) %>%
    filter(bin <= flank | bin >= (total_bins - flank)) %>%
    summarise(background = mean(signal))
}

calculate_tss_score <- function(tss) {
  total_bins <- max(tss$bin)
  binsize <- 2000 / total_bins
  middle <- total_bins / 2
  flank <- 50 / binsize
  
  tss %>%
    group_by(sample) %>%
    filter(bin >= (middle - flank) | bin <= (middle - flank)) %>%
    filter(enrichment == max(enrichment))
  
}

message("Reading deeptools profile files")
tss <- read_deeptools_profile("results/06_coverage/tss_profile.tab") %>%
  mutate(sample = str_remove(sample, "\\.raw"))

gene_body <- read_deeptools_profile("results/06_coverage/genebody_profile.tab") %>%
  mutate(sample = str_remove(sample, "\\.cpm"))

message("Calculating TSS enrichment scores")

background <- estimate_background(tss)

tss <- merge(tss, background, by = "sample") %>%
  mutate(enrichment = signal / background)

scores <- calculate_tss_score(tss)

message("Plotting TSS enrichment")
p1 <- ggplot(tss, aes(bin, enrichment)) +
  geom_vline(aes(xintercept = bin), data = scores, linetype = "dashed", color = "grey") +
  geom_line(color = "navy") +
  geom_text(aes(label = round(enrichment, 2), y = enrichment + 0.3), data = scores, size = 3) +
  ylab("TSS enrichment") +
  facet_wrap(~sample) +
  theme_bw() +
  theme(strip.background = element_blank(), panel.grid = element_blank()) +
  theme(panel.spacing.x = unit(0.5, "cm")) +
  theme(axis.title.x = element_blank(), axis.text = element_text(color = "black")) +
  scale_x_continuous(breaks = c(0, 200, 400), labels = c("-2 kb", "TSS", "2 kb"))

message("Saving to results/07_tss_enrichment/tss_enrichment.png")
ggsave("results/07_tss_enrichment/tss_enrichment.png", p1, height = 8, width = 16, units = "cm", dpi = 1200)

message("Plotting gene body coverage (CPM)")
p2 <- ggplot(gene_body, aes(bin, signal)) +
  geom_line(color = "navy") +
  ylab("CPM") +
  facet_wrap(~sample) +
  theme_bw() +
  theme(strip.background = element_blank(), panel.grid = element_blank()) +
  theme(panel.spacing.x = unit(0.5, "cm")) +
  theme(axis.title.x = element_blank(), axis.text = element_text(color = "black")) +
  scale_x_continuous(breaks = c(0, 200, 400, 600), labels = c("-2 kb", "TSS", "TES", "2 kb"))

message("Saving to results/07_tss_enrichment/gene_body.png")
ggsave("results/07_tss_enrichment/gene_body.png", p2, height = 8, width = 16, units = "cm", dpi = 1200)
