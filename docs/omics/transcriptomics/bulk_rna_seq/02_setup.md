## Setting Up The Analysis Directory


```sh
mkdir rna_seq_pipeline
cd rna_seq_pipeline
```

```sh
mkdir data
mkdir tools
mkdir qc_output
mkdir trimmed_output
mkdir star_index
mkdir alignment_output
mkdir featurecounts_output
```


## Creating A Conda Environment

```sh
cd tools
wget https://raw.githubusercontent.com/BioNomad/omicsTrain/main/docs/omics/transcriptomics/bulk_rna_seq/data/rnaseq_environment.yml
```

```sh
conda env create -f rnaseq_environment.yml    # create conda environment
source activate rnaseq                        # activate conda environment
cd ..                                         # leave tools directory
```

## Downloading Fastq Read Data

```sh
cd data
nano accList.txt
```

!!! info "accList.txt"

    ```sh
    SRR8440545
    SRR8440550
    SRR8440537
    SRR8440481
    ```

```sh
nano meta.txt
```

!!! info "meta.txt"

    ```sh
    ID  Diagnosis Age Sex
    SRR8440545  Control 53  male
    SRR8440550  Control 81  male
    SRR8440537  AD  74  female
    SRR8440481  AD  79  male
    ```
    
```sh
fastq-dump -X 9999999999999  --gzip $(<./accList.txt)
```

## Download Reference Data

!!! info "[Homo sapiens chromosome 17, GRCh38.p14 Primary Assembly](https://www.ncbi.nlm.nih.gov/nuccore/NC_000017.11)"


!!! info "Downloading the Reference FASTA File"

    ![](images/fasta_download.png)
    

!!! info "Downloading the Reference GFF3 File"

    ![](images/annotation_download.png)
    
    
    
