
library(DESeq2)
library(dplyr)
library(ggplot2)

# Params
samples_sheet <- "/data/classes/seg_epig2026/data/atac_session/00_full_data/samples.csv"
counts_file <- "/data/classes/seg_epig2026/data/atac_session/00_full_data/all_cutsites.counts.csv"
ann_file <- "resources/all_cutsites.ann.txt"

results_dir <- "results/11_differential"
dir.create(results_dir, showWarnings = F)

# read coldata
message("Reading coldata")
coldata <- read.csv(samples_sheet, row.names = "sample_title")
coldata$condition <- factor(coldata$condition)


# read counts
message("Reading counts")
counts <- read.csv(counts_file)
rownames(counts) <- paste0("all_cutsites_peak_", c(1:nrow(counts)))
counts <- counts[, rownames(coldata)]

# read annotation
message("Reading annotation")
ann <- read.delim(ann_file)

message("Creating DESeqDataSet")
# object
dds <- DESeqDataSetFromMatrix(
  countData = counts, colData = coldata, design = ~condition
)

# save object
message(paste0("Saving DDS to ", file.path(results_dir, "dds.RDS")))
saveRDS(dds, file = file.path(results_dir, "dds.RDS"))

# VST normalization
vsd <- vst(dds, blind = T)

# PCA
message("Plotting PCA")
pca <- plotPCA(vsd, intgroup = "condition", n = Inf) +
  theme_classic() +
  theme(legend.title = element_blank())

message(paste0("Saving PCA to ", file.path(results_dir, "pca.png")))
ggsave(file.path(results_dir, "pca.png"), pca, height = 8, width = 12, units = "cm", dpi = 1200)


# Differential analysis
contrast <- c("SIM6_vs_CIM7")
control <- c("CIM7")
pval <- 0.05
fc <- 1.5

dds$condition <- relevel(dds$condition, ref = control)

dds <- DESeq(dds)

res <- results(dds, name = paste0("condition_", contrast), alpha = pval)
lfc <- lfcShrink(dds = dds, res = res, coef = paste0("condition_", contrast))

# Save result
df <- lfc %>%
  as.data.frame() %>%
  mutate(peak = rownames(lfc)) %>%
  arrange(padj) %>%
  mutate(dar = case_when(
    padj < pval & log2FoldChange >= log2(fc) ~ "up",
    padj < pval & log2FoldChange <= -log2(fc) ~ "down",
    TRUE ~ "NS"
  )) %>%
  mutate(dar = factor(dar, levels = c("up", "down", "NS"))) %>%
  left_join(ann, by = join_by(peak == V4))

dars <- df %>% count(dar)
print(dars)

outfile <- file.path(results_dir, paste0(contrast, ".results.csv"))
message("Saving DESeq results to ", outfile)
write.csv(df, outfile, quote = F)

# Make volcano plot
message("Plotting volcano")
volcano <- ggplot(df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = dar), size = 0.75, alpha = 0.5) +
  scale_color_manual(values = c(2, 4, "grey")) +
  geom_text(aes(label = n, x = -2.5, y = 60), data = dars[dars$dar == "down", ], color = 4) +
  geom_text(aes(label = n, x = 2.5, y = 60), data = dars[dars$dar == "up", ], color = 2) +
  labs(x = "log2 SIM6 / CIM7", y = "-log10 Padj") +
  theme_classic() +
  theme(legend.title = element_blank())

outfile <- file.path(results_dir, paste0(contrast, ".volcano.png"))
message(paste0("Saving volcano to ", outfile))
ggsave(outfile, volcano, height = 8, width = 12, units = "cm", dpi = 1200)

# DARs by feature
message("Plotting DARs per feature")
barplot <- df %>%
  filter(dar != "NS") %>%
  count(feature2, dar) %>%
  ggplot(aes(feature2, n, fill = dar)) +
  geom_col(position = "dodge") +
  labs(x = NULL, y = "DARs") +
  scale_y_continuous(expand = c(0, 0, 0.05, 0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())

outfile <- file.path(results_dir, paste0(contrast, ".barplot.png"))
message(paste0("Saving to ", outfile))
ggsave(outfile, barplot, height = 8, width = 12, units = "cm", dpi = 1200)
