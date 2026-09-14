#!/bin/bash

# Parameters
inputs_dir="results/05_coverage"
results_dir="results/08_peak_calling"
logs_dir="logs/08_peak_calling"

samples=$(cat "config/samples.csv" | cut -d "," -f1 | tail -n +2)
conditions=$(cat "config/samples.csv" | cut -d "," -f2 | tail -n +2 | sort | uniq)

# MACS options
genome_size=118504138

echo "==== 08_peak_calling.sh ===="

mkdir -p "$results_dir" "$logs_dir"

# MACS for replicates
echo -e "\n1. MACS for replicates -----"
for sample in ${samples[@]}; do
  echo "$(date): Calling peaks for $sample"
  macs3 callpeak \
    -t "$inputs_dir/$sample.bed.gz" \
    -g "$genome_size" \
    --outdir "$results_dir" \
    -n "$sample" \
    -f BED --keep-dup all --nomodel --shift -75 --extsize 150 \
    > "$logs_dir/$sample.macs3_callpeak.log" 2>&1
done

# MACS for each condition
echo -e "\n2. MACS for conditions -----"
for condition in ${conditions[@]}; do
  echo "$(date): Calling peaks for $condition"
  macs3 callpeak \
    -t "$inputs_dir/${condition}_rep"{1,2}.bed.gz \
    -g "$genome_size" \
    --outdir "$results_dir" \
    -n "$condition" \
    -f BED --keep-dup all --nomodel --shift -75 --extsize 150 \
    > "$logs_dir/$condition.macs3_callpeak.log" 2>&1
done

# Get consensus peaks
echo -e "\n3. Consensus peaks ----"
for condition in ${conditions[@]}; do
  echo "$(date): Getting consensus peaks for $condition"
  bedtools intersect \
    -a "$results_dir/${condition}_peaks.narrowPeak" \
    -b "$results_dir/${condition}_rep"{1,2}_peaks.narrowPeak \
    -c \
    | awk -v OFS="\t" '{ if ($11 > 1) print $1, $2, $3, $4, $5, $6 }' \
    > "$results_dir/${condition}_peaks.consensus.bed" \
    2> "$logs_dir/$condition.consensus.log"
done

# Get peak union set
echo -e "\n3. Union set of peaks ----"
echo "$(date): Merging all consensus peaks"
cat "$results_dir"/*_peaks.consensus.bed \
  | sortBed \
  | mergeBed \
  | awk -v OFS="\t" '{ print $1, $2, $3, "all_peak_"NR }' \
  > "$results_dir/all_peaks.bed" \
  2> "$logs_dir/all_peaks.bedtools_merge.log"