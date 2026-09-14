#!/bin/bash

results_dir="results/03_align"
input_dir="results/02_trim/fastq"
logs_dir="logs/03_align"
threads=8
samples=$(cat "config/samples.csv" | cut -d "," -f1 | tail -n +2)

fasta="/data/classes/seg_epig2026/data/atac_session/00_raw_data/TAIR10_chr_all.fasta"

mkdir -p "$results_dir" "$logs_dir"

echo "==== 03_align.sh ===="

# 1. Index the genome

echo -e "\n1. Bwa index --------"

echo "$(date): Indexing genome"
bwa index -p "$results_dir/TAIR10" "$fasta" > "$logs_dir/TAIR10.bwa_index.log" 2>&1

# 2. Align using BWA-MEM and sorting with samtools sort
echo -e "\n2. Bwa mem --------"

for sample in ${samples[@]}; do
  echo "$(date): Aligning $sample"
  (bwa mem \
    -t $threads \
    -R "@RG\tID:$sample\tSM:$sample" \
    "$results_dir/TAIR10" \
    "$input_dir/${sample}_1.trimmed.fq.gz" \
    "$input_dir/${sample}_2.trimmed.fq.gz" \
    | samtools sort -o "$results_dir/$sample.sorted.bam" -
  ) > "$logs_dir/$sample.bwa_mem.log" 2>&1
done

# 3. Picard MarkDuplicates
echo -e "\n3. Picard MarkDuplicates --------"

for sample in ${samples}; do
  echo "$(date): Marking duplicates in $sample"
  (picard MarkDuplicates \
    -I "$results_dir/$sample.sorted.bam" \
    -M "$results_dir/$sample.markdup.txt" \
    -O "$results_dir/$sample.markdup.bam" \
    && rm "$results_dir/$sample.sorted.bam" \
  ) > "$logs_dir/$sample.picard_markduplicates.log" 2>&1
done

# 3. Samtools stats
echo -e "\n4. Samtools stats ----"

for sample in ${samples[@]}; do
  echo "$(date): Running samtools stats on $sample.markdup"
  samtools stats \
    -@ "$threads" \
    "$results_dir/$sample.markdup.bam" \
    > "$results_dir/$sample.markdup.stats.txt" \
    2> "$logs_dir/$sample.markdup.samtools_stats.log"
done

