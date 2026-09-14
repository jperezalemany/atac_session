#!/bin/bash

raw_dir="/data/classes/seg_epig2026/data/atac_session/00_raw_data/"
results_dir="results/02_trim"
logs_dir="logs/02_trim"
threads=8
samples=$(cat "config/samples.csv" | cut -d "," -f1 | tail -n +2)

echo -e "==== 02_trim.sh ===="

mkdir -p $results_dir $logs_dir

# 1. Cutadapt
echo -e "\n1. Cutadapt --------"

adapter="CTGTCTCTTATACACATCT"

mkdir -p "$results_dir/fastq" "$logs_dir"

for sample in ${samples[@]}; do
  echo "$(date): Running cutadapt over $sample"
  cutadapt \
    -j "$threads" \
    -a "$adapter" \
    -A "$adapter" \
    -q 20,20 -m 20 \
    -o "$results_dir/fastq/${sample}_1.trimmed.fq.gz" \
    -p "$results_dir/fastq/${sample}_2.trimmed.fq.gz" \
    "$raw_dir/${sample}_1.fq.gz" \
    "$raw_dir/${sample}_2.fq.gz" \
    > "$logs_dir/$sample.cutadapt.log" 2>&1
done

# 2. FastQC

echo -e "\n2. FastQC --------"

mkdir -p "$results_dir/fastqc_trimmed"

ls $results_dir/fastq/*.trimmed.fq.gz | while read fastq; do
  name=$(basename -s .trimmed.fq.gz $fastq)
  echo "$(date): Running FastQC over $name..."
  fastqc \
    --outdir "$results_dir/fastqc_trimmed" \
    "$fastq" \
    > "$logs_dir/$name.fastqc_trimmed.log" 2>&1
done

# 3. MultiQC

echo -e "\n3. MultiQC -------"

echo "$(date): Running MultiQC for $results_dir/fastqc_trimmed"

multiqc \
  --force \
  -n "$results_dir/multiqc_trimmed" \
  "$results_dir/fastqc_trimmed" \
  > "$logs_dir/multiqc_trimmed.log" 2>&1
