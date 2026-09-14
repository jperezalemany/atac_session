#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

message("\n==== 05_library_qc.R ====\n")

results_dir <- "results/05_library_qc"
dir.create(results_dir, showWarnings = F)

parse_samtools_stats <- function(file) {
  
  x <- readLines(file)
  
  get_sn <- function(name) {
    line <- x[str_detect(x, paste0("^SN\\t", fixed(name), ":"))]
    
    if (length(line) == 0)
      return(NA_real_)
    
    as.numeric(
      str_extract(
        str_remove(line, paste0("^SN\\t+", fixed(name), ":\\t*")),
        "^\\d+"
      )
    )
  }
  
  filename <- basename(file)
  
  data.frame(
    sample = str_remove(
      filename,
      "\\.(markdup|filter)\\.stats\\.txt$"
    ),
    total = get_sn("raw total sequences"),
    mapped = get_sn("reads mapped"),
    duplicated = get_sn("reads duplicated"),
    stringsAsFactors = FALSE
  )
}

parse_blacklist <- function(file) {
  
  filename <- basename(file)
  
  data.frame(
    sample = str_remove(
      filename,
      "\\.blacklist\\.txt$"
    ),
    blacklist = as.numeric(readLines(file)),
    stringsAsFactors = FALSE
  )
}

parse_insert_size <- function(file) {

    x <- readLines(file)

    # Keep only insert-size lines
    lines <- x[str_detect(x, "^IS\\t")]

    filename <- basename(file)

    sample <- str_remove(
        filename,
        "\\.filter\\.stats\\.txt$"
    )

    data.frame(
        sample = sample,
        insert_size = as.numeric(str_split_fixed(lines, "\\t", 6)[, 2]),
        count = as.numeric(str_split_fixed(lines, "\\t", 6)[, 3]),
        stringsAsFactors = FALSE
    )
}

blacklist_files <- list.files(
  "results/04_filter/",
  pattern = "\\.blacklist\\.txt$",
  full.names = TRUE
)

blacklist <- bind_rows(
  lapply(blacklist_files, parse_blacklist)
)

markdup_files <- list.files(
  "results/03_align/",
  pattern = ".markdup.stats.txt",
  full.names = TRUE
)

filter_files <- list.files(
  "results/04_filter/",
  pattern = ".filter.stats.txt",
  full.names = TRUE
)

message("Parsing samtools stats files")

markdup_stats <- bind_rows(lapply(markdup_files, parse_samtools_stats))

clean_stats <- bind_rows(lapply(filter_files, parse_samtools_stats)) %>%
  select(sample, total) %>%
  rename(clean = total)

blacklist <- bind_rows(lapply(blacklist_files, parse_blacklist))

insert_sizes <- bind_rows(lapply(filter_files, parse_insert_size)) %>%
  group_by(sample) %>%
  mutate(fraction = count / sum(count)) %>%
  ungroup()

stats <- markdup_stats %>%
  left_join(blacklist) %>%
  left_join(clean_stats) %>%
  mutate(pct_mapped = mapped / total * 100) %>%
  mutate(pct_duplicated = duplicated / mapped * 100) %>%
  mutate(pct_blacklist = blacklist / mapped * 100) %>%
  mutate(pct_clean = clean / total * 100)


# Clean reads percent
message("Plotting % of uniquely mapped reads")
p1 <- stats %>%
  mutate(discarded = total - clean) %>%
  pivot_longer(c(clean, discarded), names_to = "cat", values_to = "n") %>%
  mutate(cat = factor(cat, levels = c("discarded", "clean"))) %>%
  ggplot(aes(n/1E6, sample)) +
  geom_col(aes(fill = cat)) +
  geom_text(aes(label = paste0(round(pct_clean), "%"), x = clean/1E6), size = 3, data = stats, position = position_stack(vjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0, 0.05, 0.05)) +
  labs(x = "Reads (M)", y = NULL, fill = NULL) +
  theme_classic() +
  theme(legend.position = "top")


message("Saving to results/05_library_qc/uniquely_mapping.png")
ggsave("results/05_library_qc/uniquely_mapping.png", p1, height = 8, width = 12, units = "cm", dpi = 1200)

# Percent of blacklist
message("Plotting % of reads in blacklist")
p2 <- ggplot(stats, aes(blacklist/1E6, sample)) +
  geom_col(fill = "#E69F00") +
  geom_text(aes(label = paste(round(pct_blacklist), "%")), hjust = -0.1) +
  scale_x_continuous(expand = c(0, 0, 0.05, 0.05)) +
  labs(x = "Reads in blacklist (M)", y = NULL, fill = NULL) +
  theme_classic() +
  theme(legend.position = "top")

message("Saving to results/05_library_qc/blacklist.png")
ggsave("results/05_library_qc/blacklist.png", p2, height = 6, width = 12, units = "cm", dpi = 1200)

# Percent of duplicates
message("Plotting % of duplicated reads")
p3 <- ggplot(stats, aes(duplicated/1E6, sample)) +
  geom_col(fill = "#56B4E9") +
  geom_text(aes(label = paste(round(pct_duplicated), "%")), hjust = -0.1) +
  scale_x_continuous(expand = c(0, 0, 0.05, 0.05)) +
  labs(x = "Reads duplicated (M)", y = NULL, fill = NULL) +
  theme_classic() +
  theme(legend.position = "top")
message("Saving to results/05_library_qc/duplicated.png")
ggsave("results/05_library_qc/duplicated.png", p3, height = 6, width = 12, units = "cm", dpi = 1200)

# Insert sizes
message("Plotting insert sizes")
p4 <- ggplot(insert_sizes, aes(insert_size, fraction)) +
  geom_line() +
  facet_wrap(~sample) +
  labs(x = "Fragment size (bp)", y = "Reads fraction") +
  theme_bw() +
  theme(strip.background = element_blank(), panel.grid = element_blank())

message("Saving to results/05_library_qc/fragment_sizes.png")
ggsave("results/05_library_qc/fragment_sizes.png", p4, height = 8, width = 16, units = "cm", dpi = 1200)

message("Printing stats table")
print(stats)
