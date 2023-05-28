## ChIP-Seq Analysis Pipeline


```sh
# Step 1: Quality Control with FastQC
for file in *.fastq; do
    fastqc $file
done

# Step 2: Trim Paired End Data with Trim Galore
for file in *_R1.fastq; do
    base=$(basename $file _R1.fastq)
    trim_galore --paired $base_R1.fastq $base_R2.fastq
done

# Step 3: Align to Reference Genome with Bowtie2
for file in *_val_1.fq; do
    base=$(basename $file _val_1.fq)
    bowtie2 -x reference_genome -1 $base_val_1.fq -2 $base_val_2.fq -S $base.sam
done

# Step 4: Filter BAM Files with Samtools and Sambamba
for file in *.sam; do
    base=$(basename $file .sam)
    samtools view -bS $base.sam | sambamba view -f bam -F not unmapped -S -o $base.bam /dev/stdin
done

# Step 5: Peak Calling with MACS2
for file in *.bam; do
    base=$(basename $file .bam)
    macs2 callpeak -t $base.bam -f BAM -g genome_size -n $base
done

# Step 6: Quality Assessment with Phantompeakqualtools and ChIPQC
for file in *.bam; do
    base=$(basename $file .bam)
    Rscript phantompeakqualtools.R -s $base.bam -out $base
    Rscript chipqc.R -s $base.bam -out $base
done

# Step 7: Compare Peak Calls with Bedtools
for file in *_peaks.narrowPeak; do
    base=$(basename $file _peaks.narrowPeak)
    bedtools intersect -a $base_peaks.narrowPeak -b reference_peaks.bed -wa -u > $base_filtered_peaks.bed
done

# Step 8: Combine Peak Calls with IDR and Bedtools
idr --samples $sample1_peaks.narrowPeak $sample2_peaks.narrowPeak --input-file-type narrowPeak --output-file $sample1_vs_sample2_idr_peaks.narrowPeak
bedtools merge -i $sample1_vs_sample2_idr_peaks.narrowPeak > $sample1_vs_sample2_idr_peaks_merged.bed

# Step 9: Annotate with GREAT
for file in *_idr_peaks_merged.bed; do
    base=$(basename $file _idr_peaks_merged.bed)
    great -r $file -n $base
done

```
