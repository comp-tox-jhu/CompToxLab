## FASTQ Format

The sequencing data we will be looking at is composed of reads, or short nucleotide sequences. Read information is stored in FASTQ files:

!!! info "[FASTQ Format](https://www.drive5.com/usearch/manual7/fastq_files.html)"

    ![](images/fastq_format.png)
   
   
Here we see that for each read there are four lines:

- sequence label
- nucleotide sequence
- a separator
- quality score string

## Quality Scores

Quality scores are a score representing the probability that a base was called in error. A common score used is the Phred Score:

|Phred |Quality Score	|Probability of incorrect base call	|Base call accuracy|
|-|-|-|-|
|10	|1 in 10|	90%|
|20	|1 in 100|	99%|
|30	|1 in 1000|	99.9%|
|40	|1 in 10,000|	99.99%|
|50	|1 in 100,000|	99.999%|
|60	|1 in 1,000,000|	99.9999%|

Typically, we would like the bases in our sequence to have a quality score higher than 20, meaning that there was a 1 in 100 chance the base was called in error.

## Quality Control

To check quality scores (and a few other metrics), we can use a tool called FastQC. Let's make a script:

```sh
nano initial_qc.sh
```

!!! info "inital_qc.sh"

    ```sh
    #!/bin/bash
        
    # Step 1: Perform quality control using FastQC
    mkdir qc_output/initial_qc
    fastqc -o qc_output/initial_qc/ data/*.fastq.gz
        
    # Step 2: Aggregate quality control results using MultiQC
    multiqc qc_output/initial_qc -o qc_output/initial_qc --title "Raw Data QC Report"
        
    echo "Initial Quality Control Comlete!"
    ```

## References

1. [Base-Calling of Automated Sequencer Traces Using Phred. II. Error Probabilities](https://genome.cshlp.org/content/8/3/186.full)
2. [Quality Control](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html)
