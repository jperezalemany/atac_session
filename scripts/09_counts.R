library(DiffBind)
library(rtracklayer)
library(BiocParallel)

# Params
threads <- 8
summit <- 100
peaks <- "results/08_peak_calling/all_peaks.bed"
sample_sheet <- "config/samples.csv"

results_dir <- "results/09_counts"

message("\n ==== 09_counts.R ====\n")

# Cores to use
register(MulticoreParam(workers = threads))

# Create outdir
dir.create(results_dir, showWarnings = F)

# Read samples
samples <- read.csv(sample_sheet)[, c("sample_title", "condition")]

# Read peaks
peaks <- read.table(peaks, sep="\t", header=F)
colnames(peaks) <- c("seqnames", "start", "end", "name")
peaks <- GRanges(peaks)

# Turn into DBA sample sheet
colnames(samples) <- c("SampleID", "Condition")
samples$bamReads <- paste0("results/04_filter/", samples$SampleID, ".filter.bam")
samples$Peaks <- paste0("results/08_peak_calling/", samples$SampleID, "_peaks.narrowPeak")
samples$PeakCaller <- "narrow"

# Create object
dba <- dba(sampleSheet=samples)

# Obtain counts
message("Counting reads in peaks")
counts <- dba.count(dba, peaks=peaks, score=DBA_SCORE_READS, summits=as.integer(summit), filter=0, minCount=0)
count_matrix <- as.data.frame(dba.peakset(counts, bRetrieve=TRUE))
rownames(count_matrix) <- paste0("all_peak_", c(1:nrow(count_matrix)))

message(paste0("Saving counts to ", file.path(results_dir, "all_peaks.counts.csv")))
write.csv(count_matrix, file = file.path(results_dir, "all_peaks.counts.csv"), quote=F)

# Write summits
summits <- count_matrix[, c("seqnames", "start", "end")]
summits$name <- rownames(summits)
summits$score <- "."
summits$strand <- "."

summits$start <- summits$start + 100
summits$end <- summits$start + 1

write.csv(count_matrix, file = file.path(results_dir, "all_peaks.counts.csv"), quote=F)

message(paste0("Saving summits to ", file.path(results_dir, "all_peaks.summits.bed")))
write.table(summits, file.path(results_dir, "all_peaks.summits.bed"), sep="\t", row.names = F, col.names = F, quote = F)