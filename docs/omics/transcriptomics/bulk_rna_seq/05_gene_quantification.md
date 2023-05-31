## Gene Quantification with featureCounts

Assigning which reads are expressed in which genes essentially gives each gene a "count" per sample - or how many reads actually map to that that gene. If a read maps to a gene by as little as one base pair it is counted towards that gene. However, genes can overlap with each other, meaning a read can map to more than one gene (multi-mapping). These reads are discarded as ambiguous reads given that we can't be sure of which gene they come from:

!!! info "[featureCounts Overview](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/05_counting_reads.html)"

    ![](images/featureCounts_overview.png)

## Gene Quantification

Let's create a script to convert our gff3 file to a gft file (necessary for featureCounts), remove any lines that don't have a gene or transcript id, run featureCounts on our alignment data, and perform qc with MultiQC:

!!! info "featurecounts.sh"

    ```sh
    #!/bin/bash
    
    # Step 6: Count features using featureCounts
    
    # convert gff3 to gtf
    gffread ./data/sequence.gff3 -T -o ./data/sequence.gtf
    
    # ensure there are no empty gene ids or transcript ids
    grep -v 'gene_id ""' ./data/sequence.gtf | grep -v 'transcript_id "gene' > ./data/sequence_fixed.gtf 
    
    # run featurecounts
    featureCounts \
    -a ./data/sequence_fixed.gtf \
    -o featurecounts_output/featurecounts_results.txt \
    alignment_output/*bam
    
    # run qc on featurecounts run
    mkdir ./qc_output/featurecounts_qc 
    multiqc -o ./qc_output/featurecounts_qc --title "featureCounts QC Report" ./featurecounts_output/*.summary
    ```

!!! example "Explanation of Terms"

    **gffread**:
    
    - `-T -o ./data/sequence.gtf` make a gtf file instead of a gff3 file
    
    **grep**:
    
    - `grep -v` match the opposite of this pattern
    
    **featureCounts**:
    
    - `-a ./data/sequence_fixed.gtf` path to annotation file
    - `-o featurecounts_output/featurecounts_results.txt` path to output file
    - `alignment_output/*bam` path to our alignment input data
    
    **multiqc**:
    
    - `-o ./qc_output/featurecounts_qc ` output path
    - `--title "featureCounts QC Report"` report title
    - `./featurecounts_output/*.summary` file to use to make report
    
    
## featureCounts Quality Control

Now let's examine our featureCounts QC:

!!! info "featureCounts General Statistics: percent and number of sequences assigned"

    ![](images/featureCounts_general_stats.png)
    
    
!!! info "featureCounts General Statistics: percent and number of sequences assigned, unmapped, multi-mapped, unassigned due to no features present, and unassiged due to ambiguity in mapping"

    ![](images/featureCounts_assignment_plot.png)
    
Here we note that very few reads actually mapped to genes due to multi-mapping, where one read matches multiple genes/features. Given that this is subsampled data mapped to only one chromosome, we aren't going to take any steps to fix this. However, in your own RNA-seq analysis, most reads should be assigned to features.
    
## References

1. [featureCounts Overview](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/05_counting_reads.html)
2. [featureCounts: an efficient general purpose program for assigning sequence reads to genomic features](https://academic.oup.com/bioinformatics/article/30/7/923/232889)

