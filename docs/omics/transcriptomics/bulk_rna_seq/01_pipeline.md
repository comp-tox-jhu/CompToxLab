## Bulk RNA-seq Pipeline

```sh
#!/bin/bash

# Step 1: Perform quality control using FastQC
fastqc -o qc_output/ *.fastq.gz

# Step 2: Aggregate quality control results using MultiQC
multiqc qc_output/ -o multiqc_output/

echo "Initial Quality Control Comlete!"

# Step 3: Trim data using Trim Galore
mkdir trimmed_output/
for file in *.fastq.gz; do
    trim_galore --paired --output_dir trimmed_output/ $file
done

echo "Input Data Has Been Trimmed!"

# Step 4: Index the reference genome using STAR
mkdir star_index/
star --runMode genomeGenerate \
--genomeDir star_index/ \
--genomeFastaFiles reference.fasta

echo "Genome has been indexed!"

# Step 5: Align paired Fastq files to the reference genome using STAR in a loop
mkdir alignment_output/
for file in trimmed_output/*_val_1.fq.gz; do
    base=$(basename $file _val_1.fq.gz)
    star --genomeDir star_index/ \
    --readFilesIn $file trimmed_output/$base_val_2.fq.gz \
    --outFileNamePrefix alignment_output/"$base"_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard 
done

echo "Read Alignment Complete!"

# Step 6: Count features using featureCounts
mkdir featurecounts_output/

featureCounts -a annotation.gtf \
-o featurecounts_output/counts.txt \
alignment_output/*.sam

echo "Feature Counting Complete!"

# Step 7: Remove intermediate files
rm -r qc_output/ trimmed_output/ star_index/ alignment_output/

echo "Pipeline Completed!"
```
