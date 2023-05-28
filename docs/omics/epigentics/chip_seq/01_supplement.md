## CHIP-Seq Pipeline

```sh
#!/bin/bash

# Specify input files and directories
fastq_dir="/path/to/fastq/files"
reference_genome="/path/to/reference/genome.fa"

# Create output directories
mkdir fastqc_results
mkdir trimmed_fastq
mkdir bam_files
mkdir peak_calls
mkdir quality_assessment
mkdir peak_comparison
mkdir combined_peaks
mkdir annotated_peaks

# Loop through paired-end fastq files
for file in $fastq_dir/*_R1.fastq; do
    # Perform quality control using fastqc
    fastqc -o fastqc_results $file

    # Trim data with Trim Galore
    trimmed_fastq_base=$(basename $file _R1.fastq)
    trim_galore --paired --fastqc --output_dir trimmed_fastq \
                $fastq_dir/${trimmed_fastq_base}_R1.fastq $fastq_dir/${trimmed_fastq_base}_R2.fastq

    # Align to reference genome using bowtie2
    bowtie2 -x $reference_genome -1 trimmed_fastq/${trimmed_fastq_base}_R1_val_1.fq \
            -2 trimmed_fastq/${trimmed_fastq_base}_R2_val_2.fq -S bam_files/${trimmed_fastq_base}.sam

    # Convert SAM to BAM using samtools
    samtools view -bS bam_files/${trimmed_fastq_base}.sam > bam_files/${trimmed_fastq_base}.bam

    # Sort BAM file
    samtools sort -o bam_files/${trimmed_fastq_base}_sorted.bam bam_files/${trimmed_fastq_base}.bam

    # Index BAM file
    samtools index bam_files/${trimmed_fastq_base}_sorted.bam

    # Filter BAM file using samtools and sambamba
    samtools view -b -q 30 -F 1804 bam_files/${trimmed_fastq_base}_sorted.bam > bam_files/${trimmed_fastq_base}_filtered.bam
    sambamba markdup -r -t 4 bam_files/${trimmed_fastq_base}_filtered.bam bam_files/${trimmed_fastq_base}_dedup.bam

    # Perform peak calling with MACS2
    macs2 callpeak -t bam_files/${trimmed_fastq_base}_dedup.bam -f BAMPE -g hs -n ${trimmed_fastq_base} --outdir peak_calls

    # Quality Assessment with Phantompeakqualtools
    Rscript phantompeakqualtools/run_spp.R -c=bam_files/${trimmed_fastq_base}_dedup.bam -savp -out=quality_assessment/${trimmed_fastq_base}_QC.txt

    # Quality Assessment with ChIPQC
    Rscript ChIPQC.R -s=bam_files/${trimmed_fastq_base}_dedup.bam -r=$reference_genome -o=quality_assessment/${trimmed_fastq_base}_ChIPQC

    # Compare peak calls with bedtools
    bedtools intersect -a peak_calls/${trimmed_fastq_base}_peaks.narrowPeak -b peak_calls/other_sample_peaks.narrowPeak -wa -u > peak_comparison/${trimmed_fastq_base}_shared_peaks.bed

    # Combine peak calls with IDR and bedtools
    python idr.py --samples peak_calls/sample1_peaks.narrowPeak peak_calls/sample2_peaks.narrowPeak --output combined_peaks/combined_peaks.bed

    # Annotate peaks with GREAT
    python annotate_peaks.py -i combined_peaks/combined_peaks.bed -o annotated_peaks/${trimmed_fastq_base}_annotated.bed
done
```


```sh
Step 1: Quality Control with FastQC

bash
Copy code
for file in *.fastq; do
    fastqc "$file"
done
Step 2: Trim Paired End Data with Trim Galore

bash
Copy code
for file in *_R1.fastq; do
    base=$(basename "$file" _R1.fastq)
    trim_galore --paired "$base"_R1.fastq "$base"_R2.fastq
done
Step 3: Align to Reference Genome with Bowtie2

bash
Copy code
for file in *_val_1.fq; do
    base=$(basename "$file" _val_1.fq)
    bowtie2 -x reference_genome -1 "$base"_val_1.fq -2 "$base"_val_2.fq -S "$base".sam
done
Step 4: Filter BAM Files with Samtools and Sambamba

bash
Copy code
for file in *.sam; do
    base=$(basename "$file" .sam)
    samtools view -bS "$base".sam | sambamba view -f bam -F "not unmapped" -S -o "$base".bam /dev/stdin
done
Step 5: Peak Calling with MACS2

bash
Copy code
for file in *.bam; do
    base=$(basename "$file" .bam)
    macs2 callpeak -t "$base".bam -f BAM -g genome_size -n "$base"
done
Step 6: Quality Assessment with Phantompeakqualtools and ChIPQC

bash
Copy code
for file in *.bam; do
    base=$(basename "$file" .bam)
    Rscript phantompeakqualtools.R -s "$base".bam -out "$base"
    Rscript chipqc.R -s "$base".bam -out "$base"
done
Step 7: Compare Peak Calls with Bedtools

bash
Copy code
for file in *_peaks.narrowPeak; do
    base=$(basename "$file" _peaks.narrowPeak)
    bedtools intersect -a "$base"_peaks.narrowPeak -b reference_peaks.bed -wa -u > "$base"_filtered_peaks.bed
done
Step 8: Combine Peak Calls with IDR and Bedtools

bash
Copy code
idr --samples "$sample1_peaks.narrowPeak" "$sample2_peaks.narrowPeak" --input-file-type narrowPeak --output-file "$sample1_vs_sample2_idr_peaks.narrowPeak"
bedtools merge -i "$sample1_vs_sample2_idr_peaks.narrowPeak" > "$sample1_vs_sample2_idr_peaks_merged.bed"
Step 9: Annotate with GREAT

bash
Copy code
for file in *_idr_peaks_merged.bed; do
    base=$(basename "$file" _idr_peaks_merged.bed)
    great -r "$file" -n "$base"
done

```
