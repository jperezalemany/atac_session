#!/bin/bash

# 0. Parameters
results_dir="results/04_filter"
input_dir="results/03_align"
logs_dir="logs/04_align"

blacklist="/data/classes/seg_epig2026/data/atac_session/00_raw_data/TAIR10.blacklist.bed"
genome="/data/classes/seg_epig2026/data/atac_session/00_raw_data/TAIR10_chr_all.fasta.fai"
threads=8
samples=$(cat "config/samples.csv" | cut -d "," -f1 | tail -n +2)

mkdir -p "$results_dir" "$logs_dir"

echo "==== 04_filter.sh ===="

# 1. Generate whitelist with bedtools complement
echo -e "\n1. Whitelist (inverse blacklist) ----"

echo "$(date): Generating whitelist for TAIR10"
bedtools complement \
  -i "$blacklist" \
  -g "$genome" \
  > "$results_dir/TAIR10.whitelist.bed" \
  2> "$logs_dir/TAIR10.bedtools_complement.log"

# 2. Filter alignments with samtools view
echo -e "\n2. Filter alignments ----"

for sample in ${samples[@]}; do
  echo "$(date): Filtering $sample"
  samtools view -bh \
    -@ "$threads" \
    -f 2 -F 1804 -q 30 \
    -L "$results_dir/TAIR10.whitelist.bed" \
    "$input_dir/$sample.markdup.bam" \
    > "$results_dir/$sample.filter.bam" \
    2> "$logs_dir/$sample.samtools_view.log"
done

# Samtools index
echo -e "\n2. Index alignments ----"

for sample in ${samples[@]}; do
  echo "$(date): Running samtools index on $sample"
  samtools index \
    -@ "$threads" \
    "$results_dir/$sample.filter.bam" \
    > "$logs_dir/$sample.samtools_index.log" 2>&1
done


# 3. Samtools stats
echo -e "\n3. Samtools stats ----"

for sample in ${samples[@]}; do
  echo "$(date): Running samtools stats on $sample.filter"
  samtools stats \
    -@ "$threads" \
    "$results_dir/$sample.filter.bam" \
    > "$results_dir/$sample.filter.stats.txt" \
    2> "$logs_dir/$sample.filter.samtools_stats.log"
done

# 4. Samtools view: reads in blacklist
echo -e "\n4. Numbers of reads in blacklist ----"

for sample in ${samples[@]}; do
  echo "$(date): Running samtools view on $sample.filter"
  samtools view \
    -c \
    -@ "$threads" \
    -L "$blacklist" \
    -F 2304 \
    "$input_dir/$sample.markdup.bam" \
    > "$results_dir/$sample.blacklist.txt" \
    2> "$logs_dir/$sample.blacklist.log"
done
