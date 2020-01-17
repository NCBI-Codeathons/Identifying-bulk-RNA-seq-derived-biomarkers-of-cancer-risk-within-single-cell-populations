#find marker genes for integrated data
library(dplyr)
library(Seurat)
library(Matrix)
library(data.table)
setwd("/NYGC_hack2020/NYGC_biomarkers")

#load data
tumors.integrated <- readRDS(file = "output/tumors_integrated.rds")


#Find marker genes for each cluster and export as csv
tumors.integrated.markers <- FindAllMarkers(tumors.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tumors.integrated.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) %>% 
  write.csv(file = "output/tumors_integrated/MK_genes_TUMORS_integrated.csv")