#integration (combine individual -but similar- samples)

library(dplyr)
library(Seurat)
library(Matrix)
library(DropletUtils)
library(data.table)
library(ggplot2)
library(cowplot)
setwd("/Users/morsedb/Documents/Projects/NYGC_hack2020/NYGC_biomarkers")
#tumors.integrated <- readRDS(file = "output/tumors_integrated/tumors_integrated.rds")

#read in raw sequencing data as sparse MTX matrix
T1 <- Read10X("data/sc_hg_glioma/dense_matrices/GSM2758471", gene.column = 1)
T2 <- Read10X("data/sc_hg_glioma/dense_matrices/GSM2758472", gene.column = 1)
T3 <- Read10X("data/sc_hg_glioma/dense_matrices/GSM2758473", gene.column = 1)
T4 <- Read10X("data/sc_hg_glioma/dense_matrices/GSM2758474", gene.column = 1)
T5 <- Read10X("data/sc_hg_glioma/dense_matrices/GSM2758475", gene.column = 1)
T6 <- Read10X("data/sc_hg_glioma/dense_matrices/GSM2758476", gene.column = 1)
T7 <- Read10X("data/sc_hg_glioma/dense_matrices/GSM2758477", gene.column = 1)
T8 <- Read10X("data/sc_hg_glioma/dense_matrices/GSM2758478", gene.column = 1)

#list for piping
datasets <- list(T1,T2,T3,T4,T5,T6,T7,T8)
#initialize, filter, normalize, scale
for (i in 1:length(x = datasets)) {
  datasets[[i]] <- CreateSeuratObject(counts = datasets[[i]], project = "NYGC", min.cells = 3, min.features = 200)
  datasets[[i]][["percent.mt"]] <- PercentageFeatureSet(datasets[[i]], pattern = "^MT-")
  datasets[[i]] <- subset(datasets[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
  datasets[[i]] <- NormalizeData(datasets[[i]])
  datasets[[i]] <- FindVariableFeatures(datasets[[i]],selection.method = "vst", nfeatures = 2000)
#  datasets[[i]] <- ScaleData(datasets[[i]], vars.to.regress = "nFeature_RNA")
  datasets[[i]]$new.ident <- paste0("t", i)
  Idents(datasets[[i]]) <- 'new.ident'
}

#integrate together
tumor.anchors <- FindIntegrationAnchors(object.list = datasets, dims = 1:50)
tumors.integrated <- IntegrateData(anchorset = tumor.anchors, dims = 1:50)

# Run the standard workflow for scaling, dimensional reduction, clustering, and visualization on integrated matrix
DefaultAssay(tumors.integrated) <- "integrated"
tumors.integrated <- ScaleData(tumors.integrated, verbose = FALSE)
tumors.integrated <- RunPCA(tumors.integrated, npcs = 50, verbose = FALSE)
tumors.integrated <- FindNeighbors(tumors.integrated, dims = 1:50)
tumors.integrated <- FindClusters(tumors.integrated, resolution = 0.8)
tumors.integrated <- RunUMAP(tumors.integrated, reduction = "pca", dims = 1:50)
p1 <- DimPlot(tumors.integrated, reduction = "umap", label = TRUE, 
              repel = TRUE) + NoLegend()
p2 <- DimPlot(tumors.integrated, reduction = "umap", group.by = "new.ident")
plot_grid(p1, p2)


#save integrated analysis as rds file
saveRDS(tumors.integrated, file = "output/tumors_integrated.rds")

#export scaled data as csv
integrated_data_output <- as.data.frame(as.matrix(tumors.integrated@assays$integrated@scale.data))
fwrite(integrated_data_output, row.names = TRUE, file = "output/tumors_integrated/integrated_scaled.csv")
#export meta data
write.csv(tumors.integrated@meta.data, 'output/tumors_integrated/metadata.csv')






