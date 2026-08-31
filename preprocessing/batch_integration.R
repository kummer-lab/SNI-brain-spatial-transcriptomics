# This script integrates raw data from 4 Spatial Slides
# This script requires large RAM resources
# This script was therefore executed on an HPC of Medical University Innsbruck

library(Seurat)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyverse)
library(dplyr)
#library(STutility)
library(gridExtra)
#library(ggcorrplot)
library("pheatmap")
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
#library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)  
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
#library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)

set.seed(7)

# load the data
##################################
#meta_1 = metadata("Metadata_sham_KK3.csv",sham_data_1)
#sham_data_1@meta.data = meta_1############

metadata <- function(path, experiment){
  metadata_1 = read.csv(path)
  metadata_1$X = str_split_fixed(metadata_1$X, "_", 2)[,1]
  rownames(metadata_1) = metadata_1$X
  metadata_1 = metadata_1[, 5:6]
  
  new_meta = merge(experiment@meta.data, metadata_1, by = 0)
  rownames(new_meta) = new_meta$Row.names
  ix <- match(as.character(colnames(experiment)), as.character(rownames(new_meta)))
  new_meta <- new_meta[ix,]
  return(new_meta)
}

################################################

setwd("~/storage/spatial_transcriptomics/sni_visum_brain/spaceranger/")

####### SNI 1

# load the appropriate image first
sni_image_1 = Seurat::Read10X_Image(
  image.dir = "visium_pain_sample_1/outs/spatial",
  image.name = "tissue_hires_image.png",
  filter.matrix = TRUE
)

# function to load metadata
sni_data_1 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "visium_pain_sample_1/outs/", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "SNI1", # specify name of the stored image
  image = sni_image_1,
  filter.matrix = TRUE, 
  to.upper = FALSE
)

sni_data_1@meta.data = metadata("Metadata_SNI_KK1.csv", sni_data_1)

########## SHAM 1

sham_image_1 = Seurat::Read10X_Image(
  image.dir = "visium_pain_sample_3/outs/spatial/",
  image.name = "tissue_hires_image.png",
  filter.matrix = TRUE
)


sham_data_1 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "visium_pain_sample_3/outs/", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "sham1", # specify name of the stored image
  image = sham_image_1,
  filter.matrix = TRUE, 
  to.upper = FALSE
)

sham_data_1@meta.data = metadata("Metadata_sham_KK3.csv",sham_data_1)


##### SNI 2

sni_image_2 = Seurat::Read10X_Image(
  image.dir = "visium_pain_sample_2/outs/spatial/",
  image.name = "tissue_hires_image.png",
  filter.matrix = TRUE
)

sni_data_2 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "visium_pain_sample_2/outs/", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "SNI2", # specify name of the stored image
  image = sni_image_2,
  filter.matrix = TRUE, 
  to.upper = FALSE
)

sni_data_2@meta.data = metadata("Metadata_SNI_KK2.csv",sni_data_2)

### Sham 2
sham_image_2 = Seurat::Read10X_Image(
  image.dir = "visium_pain_sample_4/outs/spatial/",
  image.name = "tissue_hires_image.png",
  filter.matrix = TRUE
)

sham_data_2 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir = "visium_pain_sample_4/outs/", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "sham2", # specify name of the stored image
  image = sham_image_2,
  filter.matrix = TRUE, 
  to.upper = FALSE
)


sham_data_2@meta.data = metadata("Metadata_sham_KK4.csv",sham_data_2)

# add the condition and the batch effect
####################################################################

sni_data_1$condition <- "SNI"
sni_data_1$batch <- "Batch1"
sham_data_1$condition <- "sham"
sham_data_1$batch <- "Batch2"
sni_data_2$condition <- "SNI"
sni_data_2$batch <- "Batch3"
sham_data_2$condition <- "sham"
sham_data_2$batch <- "Batch4"


#####################################################################
sni_data_1 <- subset(sni_data_1, nCount_Spatial>0)
sham_data_1 <- subset(sham_data_1, nCount_Spatial>0)
sni_data_2 <- subset(sni_data_2, nCount_Spatial>0)
sham_data_2 <- subset(sham_data_2, nCount_Spatial>0)

experiment.list <- c(sni_data_1, sham_data_1, sni_data_2, sham_data_2)
for (i in 1:length(experiment.list)) {
  experiment.list[[i]] <- SCTransform(experiment.list[[i]], assay = "Spatial", vst.flavor="v2", verbose = FALSE)
}

# select the feeatures for downstream integration 
experiment.features <- SelectIntegrationFeatures(object.list = experiment.list, nfeatures = 3000)
experiment <- PrepSCTIntegration(object.list = experiment.list, anchor.features = experiment.features, 
                                 verbose = TRUE)

experiment.anchors <- FindIntegrationAnchors(object.list = experiment, normalization.method = "SCT", 
                                             anchor.features = experiment.features, verbose = TRUE, dims = 1:30)
experiment.integrated <- IntegrateData(anchorset = experiment.anchors, normalization.method = "SCT", 
                                       verbose = TRUE, dims = 1:30)



########################################################################
# label transfer
########################################################################
DefaultAssay(experiment) <- "Spatial"

experiment = NormalizeData(experiment) %>% 
  ScaleData()
markers = FindAllMarkers(experiment)
Idents(experiment) = experiment@meta.data$label
markers_per_region = FindAllMarkers(experiment)


#save.image(file = "../data_analysis/dz_Visum_Pain.RData")
save.image(file = "../data_analysis/Visum_Pain.RData")


##### Pre-Analysis checks 

# set the assay to an integrated analysis 
DefaultAssay(experiment.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
experiment <- RunPCA(experiment.integrated, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
experiment <- RunUMAP(experiment, reduction = "pca", dims = 1:25)
experiment <- RunTSNE(experiment, reduction = "pca", dims = 1:25)
experiment <- FindNeighbors(experiment, reduction = "pca", dims = 1:25)
experiment <- FindClusters(experiment)#, resolution = 0.7)


#######################################################################
p1 <- DimPlot(experiment, reduction = "umap", group.by = "condition")
p2 <- DimPlot(experiment, reduction = "umap") 
p3 <- DimPlot(experiment, reduction = "umap",group.by = "batch")
p4 <- DimPlot(experiment, reduction = "umap", group.by = "labels", label = T) + NoLegend()
SpatialDimPlot(experiment, group.by = "labels", label = T, ncol = 2, label.size = 3) + NoLegend()
SpatialDimPlot(experiment, group.by = "labels", images = c("sham1", "sham2"), label = T, ncol = 2, label.size = 3) + NoLegend()

clusters <- DimPlot(experiment, reduction = "umap", group.by = "seurat_clusters", label = T) + NoLegend()


a1 = DimPlot(experiment, reduction = "tsne", group.by = "seurat_clusters", label = T) + NoLegend()
a2 = DimPlot(experiment, reduction = "umap", group.by = "seurat_clusters", label = T) + NoLegend()
a3 = DimPlot(experiment, reduction = "tsne", group.by = "labels", label = T) + NoLegend()
a4 = DimPlot(experiment, reduction = "umap", group.by = "labels", label = T) + NoLegend()
grid.arrange(a1,a2,a3,a4,nrow = 2, ncol = 2)


DefaultAssay(experiment) = "SCT"
Idents(experiment) = "seurat_clusters"
experiment = PrepSCTFindMarkers(experiment)
seurat_cluster_degs = FindAllMarkers(experiment)

modified_degs = seurat_cluster_degs[seurat_cluster_degs$avg_log2FC>0,]
modified_degs = modified_degs %>% group_by(cluster) %>% slice(1:2) # order(p_val_adj) %>% 

VlnPlot(experiment, "Fos")
FeaturePlot(experiment, "Pvalb", label = T, split.by = "batch")

labels <- DimPlot(experiment, reduction = "umap", group.by = "labels", label = T) + NoLegend()
clusters + labels

SpatialFeaturePlot(experiment, features = c("nCount_Spatial", "nFeature_Spatial"), ncol = 4)
#p2 +p3 + p1 + p4
grid.arrange(p2, p3,p1, p4,nrow = 2, ncol = 2)

# make a barplot of cells and counts
condition_counts = experiment@meta.data %>% as.data.frame.table
counts = as.data.frame(table(condition_counts$Freq.seurat_clusters,condition_counts$Freq.condition))

final_counts = counts %>% group_by(Var1, Var2) %>% summarise(percentage = sum(Freq)) %>%
  mutate(freq = percentage/ sum(percentage))

ggplot(final_counts, aes(fill = Var2, y=freq, x=Var1, label = round(freq, 2))) +  
  geom_bar(stat="identity") +  
  geom_text(size = 3, position = position_stack(vjust = 0.5))









