import sys, getopt, subprocess
from src.common.load_h5 import H5COUNTS
from src.preprocess.build_h5_GSE103224 import build_h5
import pandas as pd
import warnings

# warnings.simplefilter("ignore")

# Load data
scRNAdata = H5COUNTS('data/GSE103224.h5')
# Preprocess data
scRNAdata.preprocess_data(log_normalize=True, filter_genes=False, n_neighbors=False, umap=False)
# # Add clustering results
scRNAdata.add_clustering_results("data/clustering_aligned/metadata.csv")


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
        opts, args = getopt.getopt(argv, "h:g:", ["ifile1=", "ifile2="])
        print(args, opts)

    except getopt.GetoptError:
        print('run.py -g <genes> -r <resolution>')
        sys.exit(3)


    for opt, arg in opts:
        if opt == '-h':
            print('python get_de_genes.py -g <genes>')
            sys.exit()
        elif opt in ("-g", "--genes"):
            biomarker = arg

    scRNAdata.get_aggregated_cluster_expression(biomarker, groupby=["cluster"], quantile_threshold=0.75)
    cluster_id = scRNAdata.cluster_exp_quantile[scRNAdata.cluster_exp_quantile > 0.0].index

    de_genes = scRNAdata.get_de_genes_by_cluster(int(cluster_id[0]))
    print(de_genes)
    return de_genes


if __name__== "__main__":
    main(sys.argv[1:])