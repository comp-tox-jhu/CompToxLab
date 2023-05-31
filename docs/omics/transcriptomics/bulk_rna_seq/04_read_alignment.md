## Read Alignment with STAR

Now that our reads are trimmed and we are left with high quality biological sequences we can try to map them back to the genome of the organism we sampled. A popular choice of aligner is STAR or Spliced Transcripts Alignment to a Reference. STAR works by finding something called a Maximal Mappable Prefix (MMP) or the longest part of a read that matches a reference sequence. STAR is splice aware in that when it cannot be mapped contiguously to a genome (like in that case where the MMP maps to an exon and and intron!) the next MMP will be identified and will be mapped to only an acceptor splice site! This MMP, when failing to reach the end of a read, will be extended till it can find regions that do match - allowing for mismatches and indels:


!!! info "STAR: ultrafast universal RNA-seq aligner"

    ![](images/star_aligner_overview.jpeg)
    
## Building an Index

Before we can map our reads we need to index our genome, in order to determine order sequence information in such a way to speed up alignment:

```sh
nano index.sh
```

!!! info "index.sh"

    ```sh
    #!/bin/bash
    
    # Step 4: Index the reference genome using STAR
    star --runMode genomeGenerate \
    --genomeDir star_index/ \
    --genomeFastaFiles data/sequence.fasta
    ```

## Gene Annotation

The choice of annotation is very important when aligning sequencing! Different gene annotation sources will have different coordinates for genes (also called features):

- **GENECODE:** 
- **UCSC or ENSEMBL:**


## Read Alignment

!!! info "align.sh"

    ```sh
    #!/bin/bash
    
    # Step 5: Align Fastq Files to a Reference Genome
    
    # unzip fastq files
    gzip -d ./trimmed_output/*.fq.gz
    
    for file in trimmed_output/*.fq; 
    do
        base=$(basename $file .fq)
        STAR --genomeDir star_index/ \
        --readFilesIn $file \
        --outFileNamePrefix alignment_output/${base}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --sjdbGTFfile ./data/sequence.gff3 \
        --sjdbOverhang 49
    done
    ```
   
## References

1. [STAR: ultrafast universal RNA-seq aligner](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/)
2. [GENECODE](https://www.gencodegenes.org/)
3. [ENSEMBL](http://www.ensembl.org/index.html?redirect=no)
4. [UCSC](https://genome.ucsc.edu/cgi-bin/hgGateway)
