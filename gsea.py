import sys, getopt, subprocess
from src.common.load_h5 import H5COUNTS
from src.preprocess.build_h5_GSE103224 import build_h5
import pandas as pd
import warnings
from src.analysis.gsea_analysis import GSEA_Analysis
from src.visualization import heatmap
# warnings.simplefilter("ignore")

# Load data
scRNAdata = H5COUNTS('data/GSE103224.h5')
# Preprocess data
scRNAdata.preprocess_data(log_normalize=True, filter_genes=False, n_neighbors=False, umap=False)
# # Add clustering results
scRNAdata.add_clustering_results("data/interim/", tumor_ids=[1, 2, 3, 4, 5, 6, 7, 8])

gsea = GSEA_Analysis(scRNAdata, threshold=0.05,)
# # Aggregate all cell expressions to find clusters with the biomarkers expressed
# scRNAdata.get_aggregated_cluster_expression(biomarkers, quantile_threshold=0.75,)
#
# # Run GSEA on all the DE genes for each cluster
# from src.analysis.gsea_analysis import GSEA_Analysis
# gsea = GSEA_Analysis(scRNAdata, path='data/interim/', threshold=0.05,) # path leads the file with the DE genes list for each cluster
# gsea.get_gsea_result()
#
# # Get the GSEA results of only the clusters which have a query biomarker expressed
# query_biomarker = ["CDC6"]
# result = gsea.get_gsea_result_by_cluster(scRNAdata.get_clusters_with_biomarker_expression(query_biomarker))
#
# # Visualize
# from src.visualization import heatmap
# heatmap(result, height=1000, width=600)


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "h:g:o", ["ifile1=", "ifile2="])
        print(args, opts)

    except getopt.GetoptError:
        print('run.py -g <genes> -r <resolution>')
        sys.exit(3)


    for opt, arg in opts:
        if opt == '-h':
            print('python gsea.py -g <genes, e.g. SLC16A3,FZD3> -o <output_file_path>')
            sys.exit()
        elif opt in ("-g", "--genes"):
            biomarker = arg

        elif opt in ("-o", "--output"):
            output = arg

    scRNAdata.get_aggregated_cluster_expression(biomarker if type(biomarker) == list else [biomarker], )
    result = gsea.get_gsea_result_by_cluster(scRNAdata.get_clusters_with_biomarker_expression(biomarker))
    heatmap(result, file_output=output)
    print("File output to", output)

if __name__== "__main__":
    main(sys.argv[1:])