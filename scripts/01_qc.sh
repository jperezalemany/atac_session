#!/bin/bash

# 0. Parameters
raw_dir="/data/classes/seg_epig2026/data/atac_session/00_raw_data/"
results_dir="results/01_qc"
logs_dir="logs/01_qc"

mkdir -p "$results_dir" "$logs_dir"

echo -e "==== 01_qc.sh ===="

# 1. FastQC

echo -e "\n1. FastQC --------"

mkdir -p "$results_dir/fastqc_raw"

ls $raw_dir/*.fq.gz | while read fastq; do
  name=$(basename -s .fq.gz $fastq)
  echo "$(date): Running FastQC over $name..."
  fastqc --outdir "$results_dir/fastqc_raw" "$fastq" > "$logs_dir/$name.fastqc_raw.log" 2>&1
done

# 2. MultiQC

echo -e "\n2. MultiQC -------"

echo "$(date): Running MultiQC for $results_dir/fastqc_raw ..."

multiqc \
  --force \
  -n "$results_dir/multiqc_raw" \
  "$results_dir/fastqc_raw" \
  > "$logs_dir/multiqc_raw.log" 2>&1
