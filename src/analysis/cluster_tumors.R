#opening libraries
library(dplyr)
library(Seurat)

#get arguments
args = commandArgs(trailingOnly=TRUE)

DATA_DIR = args[1]
TMP_DIR = args[2]

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  USER_SET_CLUS_RES = 0.8
} else {
  USER_SET_CLUS_RES = args[3]
}

#here we are opening all of the individual runs up to MARCH152019.. these are all 20 male brains

tumor1 <- Read10X(data.dir = paste0(DATA_DIR, "/GSE103224_10X/GSM2758471"), gene.column = 1)
tumor2 <- Read10X(data.dir = paste0(DATA_DIR, "/GSE103224_10X/GSM2758472"), gene.column = 1)
tumor3 <- Read10X(data.dir = paste0(DATA_DIR, "/GSE103224_10X/GSM2758473"), gene.column = 1)
tumor4 <- Read10X(data.dir = paste0(DATA_DIR, "/GSE103224_10X/GSM2758474"), gene.column = 1)
tumor5 <- Read10X(data.dir = paste0(DATA_DIR, "/GSE103224_10X/GSM2758475"), gene.column = 1)
tumor6 <- Read10X(data.dir = paste0(DATA_DIR, "/GSE103224_10X/GSM2758476"), gene.column = 1)
tumor7 <- Read10X(data.dir = paste0(DATA_DIR, "/GSE103224_10X/GSM2758477"), gene.column = 1)
tumor8 <- Read10X(data.dir = paste0(DATA_DIR, "/GSE103224_10X/GSM2758478"), gene.column = 1)

#Create objects

# Initialize the Seurat object with the raw (non-normalized data).
T1 <- CreateSeuratObject(counts = tumor1, project = "R1", min.cells = 3,min.features = 200)
T2 <- CreateSeuratObject(counts = tumor2, project = "R2", min.cells = 3,min.features = 200)
T3 <- CreateSeuratObject(counts = tumor3, project = "R3", min.cells = 3,min.features = 200)
T4 <- CreateSeuratObject(counts = tumor4, project = "R4", min.cells = 3,min.features = 200)
T5 <- CreateSeuratObject(counts = tumor5, project = "R5", min.cells = 3,min.features = 200)
T6 <- CreateSeuratObject(counts = tumor6, project = "R6", min.cells = 3,min.features = 200)
T7 <- CreateSeuratObject(counts = tumor7, project = "R7", min.cells = 3,min.features = 200)
T8 <- CreateSeuratObject(counts = tumor8, project = "R8", min.cells = 3,min.features = 200)

datasets <- list(T1,T2,T3,T4,T5,T6,T7,T8)

#here you can run all of the functions that we have  
for (i in 1:length(x = datasets)) {
  datasets[[i]] <- subset(datasets[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000)
  datasets[[i]] <- NormalizeData(datasets[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  datasets[[i]] <- FindVariableFeatures(datasets[[i]], selection.method = "vst", nfeatures = 2000)
  datasets[[i]] <- ScaleData(datasets[[i]], vars.to.regress = "nFeature_RNA")
  datasets[[i]] <- RunPCA(datasets[[i]], features = VariableFeatures(object = datasets[[i]]), npcs = 30)
  datasets[[i]] <- FindNeighbors(datasets[[i]], dims = 1:30)
  datasets[[i]] <- FindClusters(datasets[[i]], dims = 1:30, resolution = USER_SET_CLUS_RES)
}


# find markers for every cluster compared to all remaining cells, report only the positive ones
T1_MK_clusters <- FindAllMarkers(T1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T1_MK_clusters %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)%>% 
  write.csv(.,file=paste0(TMP_DIR, "MK_genes_TUMOR1.csv"))

T2_MK_clusters <- FindAllMarkers(T2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T2_MK_clusters %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)%>% 
  write.csv(.,file=paste0(TMP_DIR, "MK_genes_TUMOR2.csv"))


T3_MK_clusters <- FindAllMarkers(T3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T3_MK_clusters %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)%>% 
  write.csv(.,file=paste0(TMP_DIR, "MK_genes_TUMOR3.csv"))


T4_MK_clusters <- FindAllMarkers(T4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T4_MK_clusters %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)%>% 
  write.csv(.,file=paste0(TMP_DIR, "MK_genes_TUMOR4.csv"))

T5_MK_clusters <- FindAllMarkers(T5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T5_MK_clusters %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)%>% 
  write.csv(.,file=paste0(TMP_DIR, "MK_genes_TUMOR5.csv"))

T6_MK_clusters <- FindAllMarkers(T6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T6_MK_clusters %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)%>% 
  write.csv(.,file=paste0(TMP_DIR, "MK_genes_TUMOR6.csv"))

T7_MK_clusters <- FindAllMarkers(T7, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T7_MK_clusters %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)%>% 
  write.csv(.,file=paste0(TMP_DIR, "MK_genes_TUMOR7.csv"))

T8_MK_clusters <- FindAllMarkers(T8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T8_MK_clusters %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)%>% 
  write.csv(.,file=paste0(TMP_DIR, "MK_genes_TUMOR8.csv"))


#here you do this for the meatada
as.data.frame(T1@meta.data)%>% 
  write.csv(.,file=paste0(TMP_DIR, "T1_META.csv"))

as.data.frame(T2@meta.data)%>% 
  write.csv(.,file=paste0(TMP_DIR, "T2_META.csv"))

as.data.frame(T3@meta.data)%>% 
  write.csv(.,file=paste0(TMP_DIR, "T3_META.csv"))

as.data.frame(T4@meta.data)%>% 
  write.csv(.,file=paste0(TMP_DIR, "T4_META.csv"))

as.data.frame(T5@meta.data)%>% 
  write.csv(.,file=paste0(TMP_DIR, "T5_META.csv"))

as.data.frame(T6@meta.data)%>% 
  write.csv(.,file=paste0(TMP_DIR, "T6_META.csv"))

as.data.frame(T7@meta.data)%>% 
  write.csv(.,file=paste0(TMP_DIR, "T7_META.csv"))

as.data.frame(T8@meta.data)%>% 
  write.csv(.,file=paste0(TMP_DIR, "T8_META.csv"))

