## Setting Up The Analysis Directory

To begin we will need a space to work, let's create a directory to house all of our input data, scripts and results:

```sh
mkdir single_cell_pipeline
cd single_cell_pipeline
```

Now let's make subfolders for our data, scripts and results:

```sh
mkdir data
mkdir tools
mkdir scripts
mkdir results
```

## Downloading Data

```sh
cd data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175814/suppl/GSE175814_RAW.tar
tar -xvf GSE175814_RAW.tar
rm GSE175814_RAW.tar
```

Moving data to it's own folder and removing the prefix:

```sh
for i in *_barcodes.tsv.gz
do 
    base=$(basename $i _barcodes.tsv.gz)
    mkdir ${base}
    mv ${base}_* ${base}/
    mv ${base}/${base}_matrix.mtx.gz ${base}/matrix.mtx.gz
    mv ${base}/${base}_barcodes.tsv.gz ${base}/barcodes.tsv.gz
    mv ${base}/${base}_features.tsv.gz ${base}/features.tsv.gz  
done
```


## Cleaning Da
