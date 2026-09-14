#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

library(org.At.tair.db)
library(clusterProfiler)
library(enrichplot)

dars_file <- "results/11_differential/SIM6_vs_CIM7.results.csv"
results_dir <- "results/12_go"

dir.create(results_dir, showWarnings = F)

message("Reading dars")
res <- read.csv(dars_file) %>% arrange(padj)
contrast <- sub(".results.csv", "", basename(dars_file))

genelist <- list()

for (d in c("up", "down")) {
    genelist[[d]] <- unique(res[res$dar == d, "geneId"])
}
sapply(genelist, length)

message("Analysing GO enrichment")
go <- compareCluster(genelist, enrichGO, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP")
go <- pairwise_termsim(go)

simp <- simplify(go)

outfile <- file.path(results_dir, paste0(contrast, ".enrichgo.csv"))
message(paste0("Saving result to ", outfile))
write.csv(as.data.frame(simp), outfile)

message("Plotting GO enrichment result")
p <- dotplot(simp, showCategory = 10) +
    scale_y_discrete(labels = function(x) {strwrap(x, 100)})

outfile <- file.path(results_dir, paste0(contrast, ".enrichgo.png"))
message(paste0("Saving dotplot to ", outfile))
ggsave(outfile, p, height = 16, width = 18, units = "cm", dpi = 1200)
