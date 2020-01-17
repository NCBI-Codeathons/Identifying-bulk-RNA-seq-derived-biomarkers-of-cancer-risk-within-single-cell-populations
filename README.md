# Comb: A pipeline for identifying bulk RNA-seq derived-biomarkers of cancer risk within single-cell populations

We aim to explore cells that express risk-conferring biomarkers, for example, from this paper: https://www.nature.com/articles/s41598-019-39273-4

Some questions we could try to answer:
- Do we see expression of risk-conferring genes in these independent scRNA-seq datasets?
- If so, are all genes expressed by one population or are they distributed between populations?
- Can we characterize these "dangerous" cell populations? (e.g. do they look like some known normal cell type)
- Do these 8 tumors differ in regards to the distribution of these "dangerous" cell populations

## Team
**Matthew	Bernstein (Lead)**\
Postdoctoral Fellow at Morgridge Institute for Research\

**Paola	Correa**\
Research Associate at HHMI Janelia Research Campus

**David	Morse**\
PhD Student at the University of Cambridge and the National Institutes of Health

**Johnny	Tran**\
Machine Learning & Bioinformatics PhD-In-Training at UT Arlington

**David	Mayhew**\
Computational Bioligist at GlaxoSmithKline in Philadelphia PA



# Readme

This readme introduces our project in the Single-Cell In the Cloud Codeathon at the NY Genome Center.

Table of Content

- [Running Our Pipeline](#Running Our Pipeline)
- [Workflow](#Workflow)

  

# Running Our Pipeline

To install our pipeline and all dependencies, run

```bash
pip install git+https://github.com/NCBI-Codeathons/Identifying-bulk-RNA-seq-derived-biomarkers-of-cancer-risk-within-single-cell-populations
```


# Workflow  
Our data consisted of single-cell RNA-seq from 8 high-grade glioma patients from  https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0567-9 paper.

### 1. Data Acquision  

1. The scRNA-seq raw data was downloaded from GSE103224. It contains 8 tumor datasets with each having a few thousand cells and ~60,000 RNA expressions. The raw data should be placed in
```
data/GSE103224_RAW.tar
```
2. Uncompress the file into the folder

```
data/GSE103224_RAW
```
3. Convert the raw data format to .mtx file and a .h5 file. For python, run the following in command line
```
$ src/preprocess/build_h5_GSE103224.py
```

### 2. Data Preprocessing  
Two workflows were used, where R processes the .mtx file and Python processes the .h5 file, using Seurat and Scanpy, respectively.

```bash
$ python run.py --preprocess
```

### 3. Data Clustering  
We performed gene/cell filtering and clustering on the filtered expression matrix using Seurat v3. The resolution used for clustering was eiather 0.8 (single tumor) and 1 (integrated data).
```
Done in Seurat
```
Clustering results per tumor are saved at

```
data/interim/
```

Clustering results for aligned-tumor single-cell expression data are saved at

~~~
data/clustering_aligned/
~~~



### 4. Find differentially-expressed genes from each cluster (for each tumor)
After obtaining the clustering results, we generate a list of differentially-expressed genes for each cluster against all other clusters. In this example, we test with gene TIMP4
```
$ pythonw get_de_genes_cluster.py -g TIMP4
```

### 5. Perform Gene Set Enrichment Analysis on the DE genes for all clusters
We want to identify whether certain cell clusters have some GO terms enriched. 

Then, by comparing across all tumors, we can see if there is a common GO term involved in at least one cluster. 

To accomplish this analysis, we build a table of Adjusted P-Value (FDR rate) by running GSEA on the DE genes of every tumor-cluster pairs.

```python
from src.analysis.gsea_analysis import GSEA_Analysis

```

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>cellular response to copper ion (GO:0071280)</th>
      <th>response to copper ion (GO:0046688)</th>
      <th>cellular response to zinc ion (GO:0071294)</th>
      <th>cellular response to cadmium ion (GO:0071276)</th>
      <th>cellular zinc ion homeostasis (GO:0006882)</th>
    </tr>
    <tr>
      <th>tumor</th>
      <th>cluster</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="5" valign="top">PJ016</th>
      <th>0</th>
      <td>0.000008</td>
      <td>0.000010</td>
      <td>0.000222</td>
      <td>0.001000</td>
      <td>0.001056</td>
    </tr>
    <tr>
      <th>1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>...</th>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th rowspan="5" valign="top">PJ048</th>
      <th>5</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>6</th>
      <td>0.008055</td>
      <td>0.011066</td>
      <td>0.007067</td>
      <td>0.015403</td>
      <td>0.016894</td>
    </tr>
    <tr>
      <th>7</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>8</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>10</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>71 rows × 5 columns</p>
### 6. Find certain clusters which have cells expressing a biomarker of interest

Given the user's biomarker of interest they want to explore, we can identify the clusters which have cells expressing this biomarker. Then, they can look at the gene sets enriched by these clusters.

### 7. Gene sets enriched for clusters with SLC16A3 expressed

![SLC16A3_GSEA](/Users/jonny/Downloads/SLC16A3_GSEA.png)
