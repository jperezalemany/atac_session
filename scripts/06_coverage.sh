#!/bin/bash

raw_dir="/data/classes/seg_epig2026/data/atac_session/00_raw_data/"
genome="$raw_dir/TAIR10_chr_all.fasta.fai"
genes="resources/Araport11_genes.subset.bed"
input_dir="results/04_filter"
results_dir="results/06_coverage"
logs_dir="logs/06_coverage"
threads=8

samples=$(cat "config/samples.csv" | cut -d "," -f1 | tail -n +2)

echo "==== 05_coverage.sh ===="

mkdir -p "$results_dir" "$logs_dir"

echo -e "\n1. Tn5 shift correction ----"
# 1. bedtools shift (Tn5 correction)
for sample in ${samples[@]}; do
    echo "$(date): Correcting Tn5 cutsite for $sample"
    (bedtools bamtobed \
        -i "$input_dir/$sample.filter.bam" \
        | bedtools shift -p 4 -m -5 -g "$genome" -i - \
        | gzip > "$results_dir/$sample.bed.gz"
    ) > "$logs_dir/$sample.bedtools_shift.log" 2>&1 
done

# 2. bedtools genomecov (tracks)
echo -e "\n2. Genome coverage -----"
five_prime_awk='{
    if ($6 == "+")
        print $1, $2, $2+1, $4, $5, $6
    else
        print $1, $3-1, $3, $4, $5, $6
}'

for sample in ${samples[@]}; do
    echo "$(date): Computing coverage for $sample"
    (zcat "$results_dir/$sample.bed.gz" \
        | awk -v OFS="\t" "$five_prime_awk" \
        | bedtools slop -b 75 -g "$genome" -i - \
        | bedtools genomecov -bg -g "$genome" -i - \
        > "$results_dir/$sample.raw.bedGraph" \
    ) > "$logs_dir/$sample.bedtools_genomecov.log" 2>&1
done

# 3. CPM normalization

echo -e "\n3. CPM normalization"

echo -e "sample\tcpm_factor" > "$results_dir/libsize_factors.txt"
for sample in ${samples[@]}; do
    factor=$(cat "$input_dir/$sample.filter.stats.txt" | grep "reads mapped:" | cut -f3 | awk '{ print 1E6 / $1 }')
    echo -e "$sample\t$factor" >> "$results_dir/libsize_factors.txt"
    echo -e "$(date): Normalizing $sample by CPM"
    awk \
        -v OFS="\t" -v factor="$factor" \
        '{ print $1, $2, $3, $4*factor}' \
        "$results_dir/$sample.raw.bedGraph" \
        > "$results_dir/$sample.cpm.bedGraph" \
        2> "$logs_dir/$sample.awk_cpm.log"
done


# 4. bedGraphToBigWig
echo -e "\n4. bedGraphToBigWig"

for bedgraph in $(ls "$results_dir/"*.bedGraph); do
    echo "$(date): Converting $(basename $bedgraph) to bigwig"
    name="$(basename -s .bedGraph $bedgraph)"
    bedGraphToBigWig \
        "$results_dir/$name.bedGraph" "$genome" "$results_dir/$name.bw" \
        > "$logs_dir/$name.bedgraphtobigwig.log" 2>&1 \
        && rm "$results_dir/$name.bedGraph"
done


# 5. deepTools TSS profile

echo -e "\n5. deepTools profiles"

echo "$(date): Computing matrix for TSS profile"
computeMatrix reference-point \
    -R "$genes" \
    -S "$results_dir"/*.raw.bw \
    -b 2000 -a 2000 \
    --missingDataAsZero \
    -p "$threads" \
    -o "$results_dir/tss_matrix.gz" \
    > "$logs_dir/tss_matrix.log" 2>&1

echo "$(date): Plotting TSS profile"
plotProfile \
    -m "$results_dir/tss_matrix.gz" \
    -o "$results_dir/tss_profile.png" \
    --outFileNameData "$results_dir/tss_profile.tab" \
    > "$logs_dir/tss_profile.log" 2>&1

# 6. deepTools TSS-TES profile
echo "$(date): Computing matrix for gene body profile"
computeMatrix scale-regions \
    -R "$genes" \
    -S "$results_dir"/*.cpm.bw \
    -b 2000 -a 2000 -m 2000 \
    --missingDataAsZero \
    -p "$threads" \
    -o "$results_dir/genebody_matrix.gz" \
    > "$logs_dir/genebody_matrix.log" 2>&1

echo "$(date): Plotting gene body"
plotProfile \
    -m "$results_dir/genebody_matrix.gz" \
    -o "$results_dir/genebody_profile.png" \
    --outFileNameData "$results_dir/genebody_profile.tab" \
    > "$logs_dir/genebody_profile.log" 2>&1