## Tumor Variant Calling

```sh
#!/bin/bash

# Step 1: Define input and output variables
fastq_dir=path/to/fastq/files
reference_genome=path/to/reference_genome.fa
output_dir=path/to/output_directory

# Step 2: Perform quality control using FastQC
fastqc -o $output_dir/fastqc_reports $fastq_dir/*.fastq

# Step 3: Trim paired-end data using Trim Galore
mkdir $output_dir/trimmed_data
for file in $fastq_dir/*_R1.fastq; do
    base=$(basename $file _R1.fastq)
    trim_galore --paired \
    --output_dir $output_dir/trimmed_data \
    $fastq_dir/${base}_R1.fastq \
    $fastq_dir/${base}_R2.fastq
done

# Step 4: Index the reference genome
bwa index $reference_genome

# Step 5: Align paired-end data with BWA
mkdir $output_dir/aligned_data
for file in $output_dir/trimmed_data/*_val_1.fq; do
    base=$(basename $file _val_1.fq)
    bwa mem $reference_genome \
    $output_dir/trimmed_data/${base}_val_1.fq \
    $output_dir/trimmed_data/${base}_val_2.fq > $output_dir/aligned_data/${base}.sam
done

# Step 6: Sort and merge paired-end files
mkdir $output_dir/sorted_data
for file in $output_dir/aligned_data/*.sam; do
    base=$(basename $file .sam)
    samtools sort -o \
    $output_dir/sorted_data/${base}.bam \
    $output_dir/aligned_data/${base}.sam
    
    samtools index $output_dir/sorted_data/${base}.bam
done

# Step 7: Mark paired-end data duplicates
mkdir $output_dir/deduplicated_data
for file in $output_dir/sorted_data/*.bam; do
    base=$(basename $file .bam)
    java -jar picard.jar MarkDuplicates \
    I=$output_dir/sorted_data/${base}.bam \
    O=$output_dir/deduplicated_data/${base}_dedup.bam \
    M=$output_dir/deduplicated_data/${base}_dedup_metrics.txt
    
    samtools index $output_dir/deduplicated_data/${base}_dedup.bam
done

# Step 8: Recalibrate bases
mkdir $output_dir/recalibrated_data
for file in $output_dir/deduplicated_data/*_dedup.bam; do
    base=$(basename $file _dedup.bam)
    gatk BaseRecalibrator \
    -R $reference_genome \
    -I $output_dir/deduplicated_data/${base}_dedup.bam \
    -O $output_dir/recalibrated_data/${base}_recal_data.table \
    --known-sites known_sites.vcf
done

# Step 9: Apply recalibration
mkdir $output_dir/final_data
for file in $output_dir/recalibrated_data/*_recal_data.table; do
    base=$(basename $file _recal_data.table)
    gatk ApplyBQSR \
    -R $reference_genome \
    -I $output_dir/deduplicated_data/${base}_dedup.bam \
    -O $output_dir/final_data/${base}_recal.bam \
    --bqsr-recal-file $output_dir/recalibrated_data/${base}_recal_data.table
done

# Step 10: Call somatic variants
mkdir $output_dir/variant_calls
for file in $output_dir/final_data/*_recal.bam; do
    base=$(basename $file _recal.bam)
    gatk Mutect2 \
    -R $reference_genome \
    -I $output_dir/final_data/${base}_recal.bam \
    -O $output_dir/variant_calls/${base}_variants.vcf
done

# Step 10: Variant annotation with Funcotator
mkdir $output_dir/annotated_variants
for file in $output_dir/variant_calls/*.vcf; do
    base=$(basename $file _variants.vcf)
    gatk Funcotator \
    -R $reference_genome \
    -V $output_dir/variant_calls/${base}_variants.vcf \
    -O $output_dir/annotated_variants/${base}_annotated.vcf
done

# Remove intermediate files
rm -rf $output_dir/fastqc_reports
rm -rf $output_dir/trimmed_data
rm -rf $output_dir/aligned_data
rm -rf $output_dir/sorted_data
rm -rf $output_dir/deduplicated_data
rm -rf $output_dir/recalibrated_data
rm -rf $output_dir/final_data
rm -rf $output_dir/variant_calls
```
