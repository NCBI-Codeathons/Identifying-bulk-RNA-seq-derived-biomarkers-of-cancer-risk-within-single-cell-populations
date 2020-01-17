import sys, getopt, subprocess

from src.common.load_h5 import H5COUNTS
from src.visualization import heatmap
from src.analysis.phate_per_tumor import phate_per_tumor
from src.analysis.gsea_analysis import GSEA_Analysis

#use Rscript to call the Seurat clustering script
def seurat_clusters_script(clusteringres):
    #call the Seurat command
    seuratcommandcall = "/usr/local/bin/Rscript src/analysis/SINGLE_CLUSTERING.R %s"%clusteringres
    subprocess.call(seuratcommandcall, shell=True)

def gsea_call(genename):

	scRNAdata = H5COUNTS('data/GSE103224.h5')
	scRNAdata.preprocess_data(log_normalize=True, filter_genes=False, n_neighbors=False, umap=False)
	scRNAdata.add_clustering_results(path='data/interim/', tumor_ids=[1, 2, 3, 4, 5, 6, 7, 8])
	# Get a list of biomarkers associated to Glioma survival
	BIOMARKER_F = "data/glioma_survival_associated_genes_Fatai.csv"
	biomarkers_df = pd.read_table(BIOMARKER_F, )
	biomarkers = pd.Index(scRNAdata.GENE_NAMES) & biomarkers_df["Gene"].unique()
	# Aggregate all cell expressions to find clusters with the biomarkers expressed
	scRNAdata.get_aggregated_cluster_expression(biomarkers, quantile_threshold=0.75,)
	# Run GSEA on all the DE genes for each cluster
	from src.analysis.gsea_analysis import GSEA_Analysis
	gsea = GSEA_Analysis(scRNAdata, path='data/interim/', threshold=0.05,) # path leads the file with the DE genes list for each cluster
	gsea.get_gsea_result()
	# Query 
	query_biomarker = genename
	result = gsea.get_gsea_result_by_cluster(scRNAdata.get_clusters_with_biomarker_expression(query_biomarker))
	heatmap(result, height=1000, width=600)
    
def main(argv):
    
    #parse command line arguments
    biomarkergene=''
    genemarkers=''
    clusteringres=''
    
    try:
        opts, args = getopt.getopt(argv, "hb:b:g:r",["ifile1=","ifile2=","ifile3="])
    except getopt.GetoptError:
        print('SeuratBiomarkerWrapper.py -b <biomarker> -g <genemarkers> -r <resolution>')
        sys.exit(3)
        
    for opt, arg in opts:
        if opt == '-h':
            print('SeuratBiomarkerWrapper.py -b <biomarker> -g <genemarkers> -r <resolution')
            sys.exit()
        elif opt in ("-b", "--biomarker"):
            biomarkergene = arg
        elif opt in ("-g", "--genemarkers"):
            genemarkers = arg
        elif opt in ("-r", "--resolution"):
            clusteringres = arg
            
    #error conditions: missing input and confliciting command line arguments
    if biomarkergene=='':
        print("Error: Missing biomarker gene input")
        sys.exit(3)
        
    #split the input of gene names
    rgenes = genemarkers.split(",")

    #Run the current Seurat clustering
    seurat_clusters_script(clusteringres)
    
    #Perform gene set enrichment analysis on the biomarker gene
    gsea_call(genename)
    
    #Run phate per tumor
    phate_per_tumor(path, biomarkergene, genemarker)
    
if __name__== "__main__":
    main(sys.argv[1:])