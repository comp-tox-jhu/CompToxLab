
## WGCNA Workflow 

```r
# --- Load Libraries -----------------------------------------------------------
library(recount3)         # manipulating recount data
library(tidyverse)        # data manipulation/plotting
library(factoextra)       # clustering
library(DGCA)             # differential correlation
library(WGCNA)            # weighted gene correlation network analysis
library(patchwork)        # plot arrangement

# --- Load Data ----------------------------------------------------------------
gtex_brain  <- readRDS("./consults/ad_concord/results/gtex_brain.rds")

counts <- assay(
  gtex_brain[,gtex_brain$gtex.smtsd=="Brain - Amygdala"],
  "counts")

meta <- data.frame(
  colData(gtex_brain[,gtex_brain$gtex.smtsd=="Brain - Amygdala"]))

features <- data.frame(
  rowRanges(gtex_brain[,gtex_brain$gtex.smtsd=="Brain - Amygdala"]))

rm(gtex_brain)
gc()

# --- Set Variables ------------------------------------------------------------

# Filtering variables
filterTypes = "central"
filterCentralType = "median"
filterCentralPercentile = 0.3

# sample outlier variables
dist_method="euclidean"
hclust_method="average"

# removing sample outlier variables
cutHeight = 50
minSize = 10

# chosing power variables
powerVector=c(c(1:10), seq(from = 12, to=20, by=2))

# module building variables
power = 5
deepSplit = 2
minClusterSize = 30
hclustMethod="ward.D2"
cutHeight=0.25

# trait correlation variables
numeric_vars = "numeric.age"

# module statistics variables
numeric_var = "numeric_age"
gene_sig_thresh = 0.2
mm_thresh = 0.8

# --- Sample Outlier Removal ---------------------------------------------------

# filter out low variance genes
filt <- filterGenes(counts,
                    filterTypes = filterTypes,
                    filterCentralType = filterCentralType,
                    filterCentralPercentile = filterCentralPercentile)

# transpose so that the samples are rows and 
# the genes are columns
filt <- t(filt)

# apply ward's clustering
hc <- hclust(d = dist(filt,method = dist_method),
                 method = hclust_method)

# sample dendrogram
hc_dendro <- fviz_dend(hc)

# cut tree
clust = WGCNA::cutreeStatic(hc,
                            cutHeight = cutHeight,
                            minSize = minSize)

# filter out outliers
keepSamples <- (clust==1)
filt <- filt[keepSamples,]
meta <- meta[keepSamples,]

# --- Soft-Thresholding Power --------------------------------------------------

# pick soft thresholding power
sft = WGCNA::pickSoftThreshold(
  filt, 
  powerVector = powerVector,
  verbose = 5)
# plot power v. fit
power_v_fit <- ggplot(sft$fitIndices,
                      aes(
                        x=sft$fitIndices[,1],
                        y=-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                        label=sft$fitIndices[,1]
                      ))+
  geom_point(size=0)+
  theme_bw()+
  labs(
    x="Soft Thresholding Power",
    y="Scale Free Topology Model Fit (R^2)"
  )+
  geom_hline(yintercept = 0.9,color="red")+
  geom_text(hjust=0, vjust=0)

# plot power v. connectivity
power_v_connectivity <- ggplot(sft$fitIndices,
                               aes(
                                 x=sft$fitIndices[,1],
                                 y=sft$fitIndices[,5],
                                 label=sft$fitIndices[,1]
                               ))+
  geom_point(size=0)+
  theme_bw()+
  labs(
    x="Soft Thresholding Power",
    y="Mean Connectivity"
  )+
  geom_text(hjust=0, vjust=0)

combined <- power_v_fit|power_v_connectivity


# --- Build Modules ------------------------------------------------------------

# calculate adjacency matrix
adjacency <- adjacency(
  filt, 
  power = power)

# calculate dissimilarity matrix
TOM <- TOMsimilarity(adjacency)
rownames(TOM) <- rownames(adjacency)
colnames(TOM) <- colnames(adjacency)
dissTOM <- 1-TOM

# Call the hierarchical clustering function
geneTree <- hclust(
  as.dist(dissTOM), 
  method = hclustMethod)

# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(
  dendro = geneTree, 
  distM = dissTOM,
  deepSplit = deepSplit, 
  pamRespectsDendro = FALSE,
  minClusterSize = minClusterSize)

# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)

# Calculate eigengenes
MEList = moduleEigengenes(filt, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(
  as.dist(MEDiss),
  method = hclustMethod)

# Plot the result
# Plot the cut line into the dendrogram
plot(METree, 
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")+
  abline(h=cutHeight, 
         col = "red")
module_eigengene_dendrogram<- grDevices::recordPlot()
plot.new()

# Call an automatic merging function
merge = mergeCloseModules(
  filt, 
  dynamicColors, 
  cutHeight = cutHeight,
  verbose = 3)

# The merged module colors
mergedColors = merge$colors

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# module dendrogram plot with dynamic and merged clusters
plotDendroAndColors(
  geneTree,
  cbind(dynamicColors, mergedColors),
  c("Dynamic Tree Cut", "Merged dynamic"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE, 
  guideHang = 0.05)

module_dendrogram<- grDevices::recordPlot()
plot.new()

# Rename to moduleColors
moduleColors <- mergedColors

# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1

# --- Trait Information --------------------------------------------------------

# Define numbers of genes and samples
nGenes <- ncol(filt)
nSamples <- nrow(filt)

# Recalculate MEs with color labels
MEs0 <- WGCNA::moduleEigengenes(
  filt,
  moduleColors)$eigengenes

MEs <- WGCNA::orderMEs(MEs0)

moduleTraitCor <- WGCNA::cor(
  MEs,
  meta %>% select(numeric_vars),
  use = "p")

moduleTraitPvalue <- WGCNA::corPvalueStudent(
  moduleTraitCor, 
  nSamples)

module_trait_df <- merge(
  moduleTraitCor,
  moduleTraitPvalue,
  by="row.names")

colnames(module_trait_df) <- c(
  paste0(numeric_vars,".correlation"),
  paste0(numeric_vars,".p.value")
)

# Will display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
dev.off()
par(mar = c(6, 8.5, 3, 3))
WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                      xLabels = colnames(meta %>% 
                                           dplyr::select(numeric_vars)),
                      yLabels = names(MEs),
                      ySymbols = names(MEs),
                      colorLabels = FALSE,
                      colors = blueWhiteRed(50),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      cex.text = .6,
                      zlim = c(-1,1),
                      main = paste("Module-trait relationships"))
trait_heatmap<- grDevices::recordPlot()
plot.new()


# --- Gene/Connectivity/Module Membership Significance -------------------------

ADJ1=abs(cor(filt,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)

GS1=as.numeric(cor(meta[numeric_var],
                   filt, 
                   use="p"))
GeneSignificance=abs(GS1)

conn_gs_plot <- do.call(
  patchwork::wrap_plots, 
  lapply(unique(moduleColors), function(module) {
    restrict1 = (moduleColors==module)
    df <- data.frame(
      connectivity=Alldegrees1$kWithin[restrict1],
      gene_significance=GeneSignificance[restrict1]
    )
    mm_gs_corr <- round(
      cor(df$connectivity,
          df$gene_significance),
      3)
    mm_gs_corr_pval <- cor.test(
      df$connectivity,
      df$gene_significance)$p.value
    
    ggplot(df,
           aes(
             x=connectivity,
             y=gene_significance
           ))+
      geom_point(color=module)+
      theme_bw()+
      geom_smooth(method = "lm",
                  color=module)+
      labs(
        x="Intramodular Connectivity",
        y="Gene Significance"
      ) +
      labs(
        title=paste0(stringr::str_to_title(module),
                     " Module"),
        subtitle = paste0(
          "correlation: ",
          as.character(mm_gs_corr),
          "\n",
          "p-value: ",
          as.character(mm_gs_corr_pval))
      )
  })) 

datME=moduleEigengenes(filt,
                       moduleColors)$eigengenes
datKME=signedKME(filt,
                 datME, 
                 outputColumnName="MM.")

conn_mm_plot <- do.call(
  patchwork::wrap_plots, 
  lapply(unique(moduleColors), function(module) {
    restrict1 = (moduleColors==module)
    df <- data.frame(
      connectivity=Alldegrees1$kWithin[restrict1],
      module_membership=(datKME[restrict1, paste("MM.", module, sep="")])^6
    )
    mm_gs_corr <- round(
      cor(df$connectivity,
          df$module_membership),
      3)
    mm_gs_corr_pval <- cor.test(
      df$connectivity,
      df$module_membership)$p.value
    
    ggplot(df,
           aes(
             x=connectivity,
             y=module_membership
           ))+
      geom_point(color=module)+
      theme_bw()+
      geom_smooth(method = "lm",
                  color=module)+
      labs(
        x="Intramodular Connectivity",
        y="Module Membership "
      ) +
      labs(
        title=paste0(stringr::str_to_title(module),
                     " Module"),
        subtitle = paste0(
          "correlation: ",
          as.character(mm_gs_corr),
          "\n",
          "p-value: ",
          as.character(mm_gs_corr_pval))
      )
  })
)

sigGenes <- lapply(
  unique(moduleColors),
  function(module){
    FilterGenes= abs(GS1)> gene_sig_thresh & abs(datKME[paste0("MM.",module)])>mm_thresh
    return(colnames(filt)[FilterGenes])
  })
names(sigGenes) <- unique(moduleColors)
```

