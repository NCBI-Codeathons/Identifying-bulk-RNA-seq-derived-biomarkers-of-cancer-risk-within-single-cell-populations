# Identifying bulk RNA-seq derived-biomarkers of cancer risk within single-cell populations

We aim to explore cells that express risk-conferring biomarkers, for example, from this paper: https://www.nature.com/articles/s41598-019-39273-4

Some questions we could try to answer:
- Do we see expression of risk-conferring genes in these independent scRNA-seq datasets?
- If so, are all genes expressed by one population or are they distributed between populations?
- Can we characterize these "dangerous" cell populations? (e.g. do they look like some known normal cell type)
- Do these 8 tumors differ in regards to the distribution of these "dangerous" cell populations

# Readme

This readme introduces our project in the Single-Cell In the Cloud Codeathon at the NY Genome Center.

Table of Content

- [Background](#Background)
- [Workflow](#Workflow)
- [Requirements](#Requirements)  

# Background


# Workflow  
Our data consisted of single-cell RNA-seq from 8 high-grade glioma patients from  https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0567-9 paper.

### 1. Data Acquision  

1. The scRNA-seq raw data was downloaded from GSE103224. It contains 8 tumor datasets with each having a few thousand cells and ~60,000 RNA expressions.
2. Convert the raw data format to .mtx file and a .h5 file.

### 2. Data Preprocessing  
Two workflows were used, where R processes the .mtx file and Python processes the .h5 file, using Seurat and Scanpy, respectively.
1. 

### 3. Data Clustering  

Three methods were implemented and compared: 

1. K-means clustering
2. Spectral clustering (nearest-neighbor and Gaussian similarity kernels)
3. Hierarchical clustering

### 3. Classifying New Data
