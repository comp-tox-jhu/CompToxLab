## Sinlge Cell RNA-Seq Data Analysis

```R
# --- Single-Cell Differential Expression Analysis Workflow --------------------
# 
# --- Load Libraries -----------------------------------------------------------
library(Seurat)               # single cell processing
library(celldex)              # single cell references
library(SingleR)              # single cell classification
library(clusterProfiler)      # enrichment analysis
library(enrichplot)           # enrichment plots
library(enrichR)              # enrichment + cell classification
library(tidyverse)            # data manipulation/plotting
library(ggfortify)            # plotting
library(patchwork)            # plotting
library(reshape2)             # manipulating data
library(ggrepel)              # plotting labels
library(gridExtra)            # plot arrangement

# --- Load Data ----------------------------------------------------------------

# read in single cell counts data
sc_counts <- readRDS(
  "./consults/wrappers/cleaned_data/sc_counts.rds"
) 

# read in single cell meta data
sc_meta <- readRDS(
  "./consults/wrappers/cleaned_data/sc_meta.rds"
)
colnames(sc_counts) <- gsub("X","",colnames(sc_counts))

# --- Set Single Cell Variables ------------------------------------------------

# set single cell object creation variables
project = "asd_organoids"
min.cells = 3
min.features = 200

# set initial qc variables
split.by = "treat"
split.plot = TRUE
nFeature_RNA_min = 200
nFeature_RNA_max = 2500
percent_mt_thresh = 5

# set highly variable features variables
selection.method = "vst"
nfeatures = 2000

# clustering variables
n_dimensions = 10
chosen_res = "1.1"

# cell marker variables
only.pos = TRUE
min.pct = 0.25
logfc.threshold = 0.25
n_markers = 2
custom_markers = c("MKI67","TOP2A","BIRC5","HES1","HOPX","TNC","NEUROG1","NHLH1","EOMES","NEUROD2","SYT4","PCP4","BCL11B","FEZF2","SATB2","DLX5","GAD2","SCGN")
hpca_label_type = "label.main"
cellmarker_filt = "Brain|Cortex|Undefined|Embryonic Prefrontal Cortex"
chosen_idents = "cellmarker"

# set differential expression variables
test.use = "MAST"
condition = "treat"

# --- Create Single Cell Object ------------------------------------------------

# create Seurat single cell object
so <- CreateSeuratObject(
  counts = sc_counts, 
  project = project, 
  meta.data = sc_meta,
  min.cells = min.cells, 
  min.features = min.features)

# remove counts/meta data not in the seurat object
rm(sc_counts,sc_meta)
gc()

# --- Initial Quality Control --------------------------------------------------

# define a column for mitochondrial content
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

# Visualize QC metrics as a violin plot
init_qc_plots <- VlnPlot(
  so, 
  features = c("nFeature_RNA", 
               "nCount_RNA",
               "percent.mt"),
  split.by = split.by,
  split.plot = split.plot,
  ncol = 3)
init_qc_plots

# filter data based on qc metrics
so <- subset(so,
             subset = nFeature_RNA > nFeature_RNA_min & nFeature_RNA < nFeature_RNA_max & percent.mt < percent_mt_thresh)

# --- Normalize the Data -------------------------------------------------------

# log normalize our data
so <- NormalizeData(so)

# --- Identify Highly Variable Features ----------------------------------------

so <- FindVariableFeatures(so,
                           selection.method = selection.method,
                           nfeatures = nfeatures)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(so), 10)

# plot variable features with labels
variable_feat_plot <- VariableFeaturePlot(so)
variable_feat_plot <- variable_feat_plot + 
  ggrepel::geom_label_repel(
    aes(label=ifelse(rownames(variable_feat_plot$data)  %in% top10,
                     rownames(variable_feat_plot$data),
                     "")
    ),
    max.overlaps=20)
variable_feat_plot


# --- Scale The Data -----------------------------------------------------------

# isolate gene names
all.genes <- rownames(so)

# scale data using gene names
so <- ScaleData(so, features = all.genes)

# --- Perform Dimension Reduction ----------------------------------------------

# run pca on data
so <- RunPCA(so,
             features = VariableFeatures(object = so))

# use the elbow/pca heatmaps plots to identify the number of pcs to use
# pca heatmap
pca_heatmap <- DimHeatmap(so, 
                          dims = 1:15,
                          cells = 500,
                          balanced = TRUE)
pca_heatmap

# elbow plot
elbow_plot <- ElbowPlot(so)
elbow_plot

# --- Cluster Cells ------------------------------------------------------------

# find nearest neighbors using KNN
so <- FindNeighbors(so, dims = 1:n_dimensions)

# cluster using the Louvain algorithm
so <- FindClusters(so, resolution = seq(.4, 1.5,by=.1))

# run umap 
so <- RunUMAP(so, dims = 1:n_dimensions)

# identify single cell resolution columns
res_cols <- names(so@meta.data)[grepl("RNA_snn_res",names(so@meta.data))]
init_umaps <- list()
for(res in res_cols){
  Idents(so) <- so[[res]]
  init_umap <- DimPlot(so, reduction = "umap")+
    labs(
      title = res
    )
  init_umaps[[res]] <- init_umap
}

# examine clusters 
n <- length(init_umaps)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(init_umaps, ncol=nCol))

# choose a resolution
Idents(so) <- so[[paste("RNA_snn_res.",chosen_res,sep = "")]]
chosen_clusters_umap <- DimPlot(so, reduction = "umap")
chosen_clusters_umap

# set cluster column
so$cluster <- so[[paste("RNA_snn_res.",chosen_res,sep = "")]]

# --- Find Cluster Markers -----------------------------------------------------

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
so.markers <- FindAllMarkers(so, 
                             only.pos = only.pos, 
                             min.pct = min.pct, 
                             logfc.threshold = logfc.threshold)
so.sig.markers <- so.markers %>%
  group_by(cluster) %>%
  slice_max(n = n_markers, 
            order_by = avg_log2FC) 

# plot marker genes
marker_dotplot <- DotPlot(so, 
                          features = unique(so.sig.markers$gene), 
                          dot.scale = 8) +
  RotatedAxis()

marker_dotplot

## custom marker dot plot
custom_marker_dotplot <- DotPlot(so, 
                                 features = custom_markers, 
                                 dot.scale = 8) +
  RotatedAxis()
custom_marker_dotplot

# ---  Assign the cell types ---------------------------------------------------

# use singler to assign cell types
hpca.se <- HumanPrimaryCellAtlasData()
so_counts <- log2(as.matrix(GetAssayData(so))+1)
pred <- SingleR(
  test = so_counts,
  ref = hpca.se, 
  assay.type.test=1,
  labels = hpca.se[[hpca_label_type]])

# add in new labels
so$singler <- pred$labels

# use cellmarker+enrichr to get cell classification
cluster_gene_list <- so.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, 
            order_by = avg_log2FC)
cluster_gene_list <- split(
  cluster_gene_list$gene,
  cluster_gene_list$cluster
)

# run through enrichr
enriched_cluster_list <- lapply(
  cluster_gene_list,
  function(x){
    enriched=enrichR::enrichr(x,databases = c(
      enrichR::listEnrichrDbs()$libraryName[grepl("CellMarker",enrichR::listEnrichrDbs()$libraryName)]
    ));
    return(enriched)
  }
)

# unlist the results and take the top result as the assigned cell_type
cellmarker <- unlist(
  lapply(
    enriched_cluster_list,
    function(x){
      unlist(x[[1]][grepl(cellmarker_filt,x[[1]]$Term),]$Term[1])}
    ))

# convert cell types to identities
names(cellmarker) <- levels(so)
so <- RenameIdents(so, cellmarker)

# add in a column for cellmarker labels and change identity back to 
# chosen resolution
so$cellmarker <- Idents(so)
Idents(so) <- so[[paste("RNA_snn_res.",chosen_res,sep = "")]]


# choose cell identities
if(chosen_idents == "singler"){
  Idents(so) <- so$singler
}else if (chosen_idents == "cellmarker"){
  Idents(so) <- so$cellmarker
}else if (chosen_idents == "cluster_numbers"){
  Idents(so) <- so[[paste("RNA_snn_res.",chosen_res,sep = "")]]
}

# plot chosen cell identities
chosen_cell_type_umap <- DimPlot(so, reduction = "umap", pt.size = 0.5) 
chosen_cell_type_umap

# --- Differential Expression --------------------------------------------------

# pairwise differential expression between individual cell types

# isolate combinations of cell types
cell_combinations <- expand.grid(levels(so),levels(so))
cell_combinations <- cell_combinations %>%
  filter(Var1!=Var2)
cell_combinations <- unique(apply(cell_combinations, 1, function(x) paste0(sort(x), collapse = "_")))

# loop through cell type combinations
# perform differential expression
between_cell_types <- list()
for (cell_comparison in cell_combinations){
  comp1 = strsplit(cell_comparison,"_")[[1]][1]
  comp2 = strsplit(cell_comparison,"_")[[1]][2]
  de <-FindMarkers(so,  
                   ident.1 = comp1,
                   ident.2 = comp2,
                   min.pct = min.pct, 
                   logfc.threshold = logfc.threshold,
                   test.use = test.use)
  between_cell_types[[cell_comparison]] = de
}

# differential expression between a cell type and everything else
de_all <-FindAllMarkers(so,  
                        min.pct = min.pct, 
                        logfc.threshold = logfc.threshold)

# pairwise differential expression between cell types in different conditions
condition_cell <- paste(unlist(so[[condition]]),unlist(so[[chosen_idents]]),sep = "-")
Idents(so) <- condition_cell
cond_cell_combinations <- expand.grid(levels(so),levels(so))
cond_cell_combinations <- cond_cell_combinations %>%
  filter(Var1 != Var2)
cond_cell_combinations <- unique(apply(cond_cell_combinations, 1, function(x) paste0(sort(x), collapse = "_")))

# loop through cell type combinations
# perform differential expression
cond_between_cell_types <- list()
for (cond_cell_comparison in cond_cell_combinations){
  comp1 = strsplit(cond_cell_comparison,"_")[[1]][1]
  comp2 = strsplit(cond_cell_comparison,"_")[[1]][2]
  de <-FindMarkers(so,  
                   ident.1 = comp1,
                   ident.2 = comp2,
                   min.pct = min.pct, 
                   logfc.threshold = logfc.threshold,
                   test.use = test.use)
  cond_between_cell_types[[cond_cell_comparison]] = de
}

# change identities back
Idents(so) <- so[[chosen_idents]]

# general bulk condition differential expression
Idents(so) <- so[[condition]]
general_bulk_de <-FindMarkers(so,  
                 ident.1 = levels(so)[1],
                 ident.2 = levels(so)[2],
                 min.pct = min.pct, 
                 logfc.threshold = logfc.threshold,
                 test.use = test.use)

# change identities back
Idents(so) <- so[[chosen_idents]]

```
