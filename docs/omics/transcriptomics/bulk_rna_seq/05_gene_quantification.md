## Gene Quantification with featureCounts


## Gene Quantification

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
