## Labelled Proteomics Analysis 

```sh
# Step 1: Convert TMT (.raw) files to mzML using OpenMS
mkdir mzml_files
for file in *.raw; do
    OpenSwathConverter -in $file -out mzml_files/${file%.raw}.mzML
done


# Step 2: Identify peptides using MSGFPlus
mkdir results
for mzml_file in mzml_files/*.mzML; do
    base_name=$(basename $mzml_file .mzML)
    MSGFPlus -s $mzml_file -o results/${base_name}_identifications.mzid -d database.fasta -t 20ppm -m 1
done
```

```R
# Step 3: Load Id files as a table
library("mzID")
mzids <- list.files(system.file('extdata', package = 'mzID'),
                    pattern = '*.mzid', full.names = TRUE)
# load in data
ids <- mzID(mzids)

# flatten into a table
fid <- flatten(id)

```
