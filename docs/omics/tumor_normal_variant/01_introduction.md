## Background

## Copy Over Data

```sh
# make a directory for the Analysis
mkdir tn_variant_analysis

# enter directory and copy over data
cd tn_variant_analysis
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz
cd ..

# make a directory for the reference
mkdir ref_data
cd ref_data
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz
```

## Initial Data QC

```sh
# make the raw data qc folder 
mkdir raw_data_qc

# run fastqc
fastqc ./raw_data/ -o raw_data_qc

# run multiqc on data 
multiqc ./raw_data_qc/ -o raw_data_qc/ --title "Raw Data QC Report"
```

## Trimmed Data 

```sh
# make an output directory
mkdir trimmed_data

# run trim galore
parallel --xapply trim_galore --paired -o ./trimmed_data/ ::: ./raw_data/*_R1.fastq.gz ::: ./raw_data/*_R2.fastq.gz
```


## Trimmed Data QC


```sh
# make the raw data qc folder 
mkdir trimmed_data_qc

# run fastqc
fastqc ./trimmed_data/*.fq.gz -o trimmed_data_qc

# run multiqc
multiqc ./trimmed_data_qc/ -o trimmed_data_qc/ --title "Trimmed Data QC Report"
```


## Index Reference

```sh
# index the reference fasta
cd ref_data
bwa index hg19.chr5_12_17.fa.gz
cd ..
```

## Read Alignment

```sh
# make alignment directory
mkdir aligned_data
cd aligned_data

for fname in ../trimmed_data/*_R1.fq.gz
do
    base=${fname%_R1*}
    bwa mem -t 4 -T 0 -R "@RG SOMETHING"../ref_data/hg19.chr5_12_17.fa.gz "../trimmed_data/${base}_R1.fq.gz"  "../trimmed_data/${base}_R2.fq.gz"  | samtools view -Shb -o ./${base}.bam 
done
```

## Sort Sam Files


```sh
cd aligned_data

# for each bam file sor the bam file by coordinate
for fbam in ./*.bam
do
    java -jar picard.jar SortSam \
    CREATE_INDEX=true \
    INPUT=fbam \
    OUTPUT=fbam.sorted.bam \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=STRICT
done
```


## Mark Duplicates

```sh
cd aligned_data

# mark duplicates
for fsortbam in ./*.sorted.bam
do
    picard -Xmx7g MarkDuplicates \
        I=./fsortbam \
        O=./fsortbam.dup.bam \
        METRICS_FILE=./fsortbam_dup_metrics.txt
done
```


## Base Recalibration

```sh

for fdupbam in ./*.dup.bam
do
    # step 1  - Build the model
    gatk --java-options "-Xmx7g" BaseRecalibrator \
    -I ./fdupbam \
    -R ../ref_data/hg19.chr5_12_17.fa.gz \
    --known-sites ../ref_data/dbsnp_146.hg38.vcf.gz \
    -O ./fdupbam_recal_data.table
    
    # step 2: Apply the model to adjust the base quality scores
    gatk --java-options "-Xmx7g" ApplyBQSR \
    -I ./fdupbam \
    -R ../ref_data/hg19.chr5_12_17.fa.gz \
    --bqsr-recal-file ./fdupbam_recal_data.table \
    -O ./fdupbam.sort.dup.bqsr.bam    
done
```

## Running the HaplotypeCaller

```sh
for fbqsr in ./*.sort.dup.bqsr.bam 
do
       
done
```
