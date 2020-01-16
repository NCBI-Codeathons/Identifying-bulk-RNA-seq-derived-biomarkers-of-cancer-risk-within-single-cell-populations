import gseapy as gp
import pandas as pd

import plotly.express as px
from ..common.load_h5 import H5COUNTS

class GSEA_Analysis():
    def __init__(self, h5_counts:H5COUNTS, path="data/interim/", threshold=0.05, tumor_ids=[1, 2, 3, 4, 5, 6, 7, 8]):
        self.gsea_table = pd.DataFrame()

        for tumor_id in tumor_ids:
            de_genes_tumor_df = pd.read_csv(path+"MK_genes_TUMOR{}.csv".format(tumor_id))
            de_genes_by_cluster = de_genes_tumor_df.groupby("cluster")["gene"].apply(lambda x: "|".join(x.unique()))

            for cluster in de_genes_by_cluster.index:
                DE_gene_list = de_genes_by_cluster[cluster].split("|")
                tumor_cluster = h5_counts.id2tumor[tumor_id] + "_" + str(cluster)
                print("Running GSEA for tumor", tumor_id)
                enr = gp.enrichr(gene_list=DE_gene_list,
                                 gene_sets=['GO_Biological_Process_2018',
                                            'GO_Cellular_Component_2018',
                                            'GO_Molecular_Function_2018'],
                                 no_plot=True,
                                 cutoff=0.05  # test dataset, use lower value from range(0,1)
                                 )
                if threshold:
                    enr.results = enr.results[enr.results["Adjusted P-value"] < threshold]
                enr_results = enr.results.set_index("Term")

                for geneset in enr_results.index:
                    self.gsea_table.loc[geneset, tumor_cluster] = enr_results.loc[geneset, "Adjusted P-value"]

    def get_gsea_result(self):
        return self.gsea_table