
## Variant Calling In Tumor Samples

```sh
#!/bin/bash

# Set input and output variables
fastq1="path/to/tumor_sample_R1.fastq.gz"
fastq2="path/to/tumor_sample_R2.fastq.gz"
ref_genome="path/to/reference_genome.fasta"
output_dir="path/to/output_directory"

# Step 1: Quality Control with FastQC
fastqc -o $output_dir $fastq1 $fastq2

# Step 2: Trimming with Trimmomatic
trimmomatic PE -threads 4 $fastq1 $fastq2 \
    $output_dir/trimmed_tumor_sample_R1.fastq.gz \
    $output_dir/unpaired_tumor_sample_R1.fastq.gz \
    $output_dir/trimmed_tumor_sample_R2.fastq.gz \
    $output_dir/unpaired_tumor_sample_R2.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 3: Indexing Reference Genome with BWA
bwa index $ref_genome

# Step 4: Alignment with BWA
bwa mem -t 4 $ref_genome \
    $output_dir/trimmed_tumor_sample_R1.fastq.gz \
    $output_dir/trimmed_tumor_sample_R2.fastq.gz \
    > $output_dir/aligned_tumor_sample.sam

# Step 5: Sorting and Merging Alignment Files with Samtools
samtools sort -@ 4 -o $output_dir/sorted_tumor_sample.bam $output_dir/aligned_tumor_sample.sam
samtools merge -@ 4 $output_dir/merged_tumor_sample.bam $output_dir/sorted_tumor_sample.bam

# Step 6: Mark Duplicates with Picard
java -jar picard.jar MarkDuplicates \
    I=$output_dir/merged_tumor_sample.bam \
    O=$output_dir/markduplicates_tumor_sample.bam \
    M=$output_dir/markduplicates_metrics.txt

# Step 7: Base Recalibration with GATK
gatk BaseRecalibrator \
    -R $ref_genome \
    -I $output_dir/markduplicates_tumor_sample.bam \
    --known-sites known_sites.vcf \
    -O $output_dir/recalibration_report.table

# Step 8: Apply Base Recalibration with GATK
gatk ApplyBQSR \
    -R $ref_genome \
    -I $output_dir/markduplicates_tumor_sample.bam \
    --bqsr-recal-file $output_dir/recalibration_report.table \
    -O $output_dir/recalibrated_tumor_sample.bam

# Step 9: Variant Calling with GATK
gatk Mutect2 \
    -R $ref_genome \
    -I $output_dir/recalibrated_tumor_sample.bam \
    -tumor tumor_sample \
    -O $output_dir/tumor_sample_variants.vcf

# Step 10: Variant Annotation with Funcotator
gatk Funcotator \
    -R $ref_genome \
    -V $output_dir/tumor_sample_variants.vcf \
    -O $output_dir/annotated_tumor_sample_variants.vcf \
    --data-sources-path funcotator_data_sources

# Step 11: Removing Intermediate Files
rm $output_dir/aligned_tumor_sample.sam
rm $output_dir/sorted_tumor_sample.bam
rm $output_dir/merged_tumor_sample.bam
rm $output_dir/markduplicates_tumor_sample.bam
rm $output_dir/recalibration_report.table

echo "Pipeline Finished Running!"

```
