## 16S Amplicon Sequencing Data Analysis

```R
# --- Amplicon Sequencing Analyis ----------------------------------------------

# --- Load Libraries -----------------------------------------------------------

# load our libraries
.libPaths(c('/cluster/tufts/hpc/tools/R/4.0.0',.libPaths()))
library(dada2)
library(phyloseq)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(phangorn)
library(msa)

# --- Quality Control ----------------------------------------------------------
# path to files
path <- "../data/raw_fastq"

# sort our files by forward and reverse strands 
# so that the sample names for each strand matches
# our data has the pattern "_pass_1.fastq.gz" 
# and "_pass_2.fastq.gz"
path2Forward <- sort(
  list.files(
    path,
    pattern="_pass_1.fastq.gz",
    full.names = TRUE)
)
path2Reverse <- sort(
  list.files(
    path,
    pattern="_pass_2.fastq.gz",
    full.names = TRUE)
)

# now let's grab our sample names
sampleNames <- sapply(
  strsplit(
    basename(path2Forward), "_"), `[`, 1)

# plot the forward strand quality plot of our first two sample
dada2::plotQualityProfile(path2Forward[1:2])+
  guides(scale = "none")

# plot the reverse strand quality plot of our first two sample
dada2::plotQualityProfile(path2Reverse[1:2])+
  guides(scale = "none")

# --- Trimming -----------------------------------------------------------------

# create new file names for filtered forward/reverse fastq files
# name each file name in the vector with the sample name
# this way we can compare the forward and reverse files 
# when we filter and trim
filtForward <- file.path(path, "filtered", paste0(sampleNames, "_F_filt.fastq.gz"))
filtReverse <- file.path(path, "filtered", paste0(sampleNames, "_R_filt.fastq.gz"))
names(filtForward) <- sampleNames
names(filtReverse) <- sampleNames

# Now we will filter and trim our sequences
out <- filterAndTrim(
  path2Forward,
  filtForward,
  path2Reverse, 
  filtReverse,
  truncLen = c(200,150),
  maxN=0, 
  maxEE=c(2,2), 
  truncQ=2, 
  rm.phix=TRUE,
  compress=TRUE)

# --- DADA2 Error Model --------------------------------------------------------
# Learn Error Rates

# dada2 uses a parametric model to learn the error rates
# for each sequence
errForward <- learnErrors(filtForward)
errReverse <- learnErrors(filtReverse)

# plot the error rate against theoretical error rates
plotErrors(errForward,nominalQ=TRUE)

# --- Inferring Sequence Variants ----------------------------------------------
# Infer Sequnce Variants

# we will now run the dada2 algorithm 
# this algorithm delivers "true" sequence variants
# with information gathered from the error model 
# generated above
dadaForward <- dada(filtForward, err=errForward)
dadaReverse <- dada(filtReverse, err=errReverse)

# let's get a summary of our first sample
dadaForward[[1]]

# --- Merging Reads ------------------------------------------------------------

# Merge Read Pairs

# so far we have "denoised", so to speak, 
# these sequence variants. We now need to merge the
# forward and reverse strands
mergers <- mergePairs(
  dadaForward,
  filtForward,
  dadaReverse, 
  filtReverse, 
  verbose=TRUE)

# --- ASV Table ----------------------------------------------------------------
# Making a Sequence Table

# now that we have merged sequences we can construct
# an Amplicon Sequence Variant (ASV) table
seqtab <- makeSequenceTable(mergers)

# --- Chimera Removal ----------------------------------------------------------

# Removing Chimeras

# Chimeric sequences occur as errors during PCR 
# when two unrelated templates for a hybrid sequence
# we will need to remove them before going forward

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE)
## check to see if the dimensions are different
## between the chimera filtered and unfiltered
## ASV tables

dim(seqtab)
dim(seqtab.nochim)

# --- Pipeline Quality Control -------------------------------------------------
# Final QC

## we have performed quite a few steps 
## and it would be nice to get a final qc check 
## before assigning taxonomy
getN <- function(x) sum(getUniques(x))
finalQC <- cbind(
  out, 
  sapply(dadaForward, getN),
  sapply(dadaReverse, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim))
colnames(finalQC) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(finalQC) <- sampleNames
finalQC

# --- Assigning Taxonomy -------------------------------------------------------
# Assigning Taxonomy

# dada2 uses a naive Bayes classifier when
# assigning taxonomy. This means we need a training
# set of sequences with known taxonomy information.
# here we use the silva database

taxa <- assignTaxonomy(seqtab.nochim, "../data/silva_nr99_v138.1_train_set.fa.gz")

# --- Constructing the Phylogenetic Tree ---------------------------------------
# extract sequences
# name the sequences with their sequence so 
# that the ends of the phylogenetic tree are labeled
# align these sequences
seqs <- getSequences(seqtab)
names(seqs) <- seqs 
mult <- msa(seqs, method="ClustalW", type="dna", order="input")

# convert multiple sequence alignment to a phyDat object
# calculate the nucleotide distances between ASVs
# use a neighbor joining algorithm to generate the tree
# finally calculate the likelihood of the tree given the sequence alignment
phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit = pml(treeNJ, data=phang.align)

# --- Making a PhyloSeq Object -------------------------------------------------

# Create phyloseq object

# upload meta data for study
# ensure the rownames of our meta data are our sample name
meta <- read.csv("../data/metaData.txt")
rownames(meta) <- meta$Run

# combine the ASV table, the meta data, and taxonomic data
# to create the phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa),
               phy_tree(fit$tree)
)

# Update ASV names to be shorter

# The full ASV DNA sequence can be hard to look at
# for this reason we move the sequence information to 
# the refseq slot of the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))


# --- Diversity Analysis -------------------------------------------------------

# Plotting Alpha Diversity Metrics
plot_richness(ps, x="Host", measures=c("Shannon", "Simpson"), color="Host")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=65,hjust=1))

# calculate the unifrac distance between samples 
# plot unifrac distances
ordu = ordinate(ps, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps, ordu, color="Host")+
  theme_bw()+
  labs(title = "Unifrac Distances")

# --- Phylum Present -----------------------------------------------------------

# transform the sample counts to proportions
# separate out our proportions
# separate our our tax info
ps.prop <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
otu = data.frame(t(data.frame(ps.prop@otu_table)))
tax = data.frame(ps.prop@tax_table) 

# merge the otu table and phylum column
# reshape our data to be accepted by ggplot
# merge taxa data with sample meta data
merged <- merge(otu,
                tax %>% select(Phylum),
                by="row.names") %>%
  select(-Row.names) %>%
  reshape2::melt() %>%
  merge(.,
        data.frame(ps.prop@sam_data) %>%
          select(Run,Host),
        by.x="variable",
        by.y="Run")

# plot our taxa 
taxa_plot <- ggplot(merged,aes(x=variable,y=value,fill=Phylum)) +
  geom_bar(stat='identity') +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  labs(
    x="",
    y="Abundance",
    title = "Barplot of Phylum Abundance"
  )+
  facet_wrap(Host ~ ., scales = "free_x")
taxa_plot

# --- Differential Abundance ---------------------------------------------------

# Differential Abundance

## convert phyloseq object to DESeq object this dataset was downsampled and 
## as such contains zeros for each ASV, we will need to
## add a pseudocount of 1 to continue and ensure the data are still integers
## run DESeq2 against Host status, and ensure wild type is control,
## filter for significant changes and add in phylogenetic info
dds = phyloseq_to_deseq2(ps, ~ Host)
dds@assays@data@listData$counts = apply((dds@assays@data@listData$counts +1),2,as.integer)
dds = DESeq(dds, test="Wald", fitType="parametric")
res = data.frame(
  results(dds,
          cooksCutoff = FALSE, 
          contrast = c("Host","C57BL/6NTac","Mus musculus domesticus")))
sigtab = res %>%
  cbind(tax_table(ps)[rownames(res), ]) %>%
  dplyr::filter(padj < 0.05) 

## order sigtab in direction of fold change
sigtab <- sigtab %>%
  mutate(Phylum = factor(as.character(Phylum), 
                         levels=names(sort(tapply(
                           sigtab$log2FoldChange, 
                           sigtab$Phylum, 
                           function(x) max(x)))))
  )

# as a reminder let's plot our abundance data again
taxa_plot

## plot differential abundance
ggplot(sigtab , aes(x=Phylum, y=log2FoldChange, color=padj)) + 
  geom_point(size=6) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle("Mus musculus domesticus v. C57BL/6NTac")

```
