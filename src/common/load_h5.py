import h5py
import numpy as np
from collections import defaultdict
import pandas as pd

from anndata import AnnData
import scanpy as sc

class H5COUNTS():
    def __init__(self, DATA_F):
        self.DATA_F = DATA_F
        with h5py.File(DATA_F, 'r') as f:
            self.CELLS = [
                str(x)[2:-1]
                for x in f['cell'][:]
            ]
            self.TUMORS = [
                str(x)[2:-1]
                for x in f['tumor'][:]
            ]
            self.GENE_IDS = [
                str(x)[2:-1]
                for x in f['gene_id'][:]
            ]
            self.GENE_NAMES = [
                str(x)[2:-1]
                for x in f['gene_name'][:]
            ]

        # Map each cell to its index in the data matrix
        self.CELL_TO_INDEX = {
            cell: index
            for index, cell in enumerate(self.CELLS)
        }

        # Map each tumor to its indices in the data matrix
        self.TUMOR_TO_INDICES = defaultdict(lambda: [])
        for index, tumor in enumerate(self.TUMORS):
            self.TUMOR_TO_INDICES[tumor].append(index)
        self.TUMOR_TO_INDICES = dict(self.TUMOR_TO_INDICES)

        self.TUMORS = np.unique(self.TUMORS)
        self.tumor2id = {tumor: tumor_id + 1 for tumor_id, tumor in enumerate(self.TUMORS)}
        self.id2tumor = {tumor_id + 1: tumor for tumor_id, tumor in enumerate(self.TUMORS)}

    def preprocess_data(self, log_normalize=True, filter_genes=False, n_neighbors=False, umap=False):
        self.tumor_to_ad = {}
        for tumor in self.TUMORS:
            print(tumor)
            counts, cells = self.counts_matrix_for_tumor(tumor)
            ad = AnnData(
                X=counts,
                obs=pd.DataFrame(data=cells, columns=['cell']),
                var=pd.DataFrame(
                    data=self.GENE_NAMES,
                    index=self.GENE_NAMES,
                    columns=['gene_name']
                )
            )
            if log_normalize:
                sc.pp.normalize_total(ad, target_sum=1e6)
                sc.pp.log1p(ad)

            if filter_genes:
                sc.pp.filter_genes(ad, min_cells=3, min_counts=200)

            if n_neighbors:
                sc.pp.neighbors(ad)

            if umap:
                sc.tl.umap(ad)

            ad.X = pd.DataFrame(ad.X, index=cells, columns=self.GENE_NAMES)
            self.tumor_to_ad[tumor] = ad

        for tumor in self.TUMORS:
            self.tumor_to_ad[tumor].X.index = pd.MultiIndex.from_tuples(self.tumor_to_ad[tumor].X.index.str.split("_", expand=True),
                                                                        names=["tumor", "cell"])
        self.all_tumor_df = pd.concat([self.tumor_to_ad[tumor].X for tumor in self.TUMORS])

    def counts_matrix_for_tumor(self, tumor):
        indices = self.TUMOR_TO_INDICES[tumor]
        with h5py.File(self.DATA_F, 'r') as f:
            counts = f['count'][indices]
        cells = list(np.array(self.CELLS)[indices])
        return counts, cells

    def get_unique_tumors(self):
        return np.unique(self.TUMORS)

    def add_clustering_results(self, path="data/interim/", tumor_ids=[1, 2, 3, 4, 5, 6, 7, 8]):
        self.tumor_cell_cluster = pd.DataFrame(columns=["tumor", "cell", "cluster"])

        for tumor in tumor_ids:
            cell_cluster = pd.read_csv(path+"T{}_META.csv".format(tumor))

            cell_cluster.rename(columns={"Unnamed: 0": "cell", "RNA_snn_res.0.8": "cluster"},
                                inplace=True)
            cell_cluster["cell"] = cell_cluster["cell"].str.split("_", expand=True)[1].astype(int) - 1
            cell_cluster["cell"] = cell_cluster["cell"].astype(str)
            cell_cluster["tumor"] = self.id2tumor[tumor]

            self.tumor_cell_cluster = self.tumor_cell_cluster.append(cell_cluster[["tumor", "cell", "cluster"]])


        self.tumor_cell_cluster["cluster"] = self.tumor_cell_cluster["cluster"].astype(str)
        self.tumor_cell_cluster.set_index(["tumor", "cell"], inplace=True)

        tumor_cell_cluster_df = pd.DataFrame(self.tumor_cell_cluster.values,
                                             index=self.tumor_cell_cluster.index.to_flat_index().str.join("_"),
                                             columns=["cluster"])
        # Assign cluster assignment to the AnnData's
        for tumor in tumor_ids:
            self.tumor_to_ad[self.id2tumor[tumor]].obs = self.tumor_to_ad[self.id2tumor[tumor]].obs.join(tumor_cell_cluster_df, on="cell")

    def get_aggregated_cluster_expression(self, biomarkers:list, quantile_threshold=0.75):
        all_tumor_cell_biomarkers = self.all_tumor_df[biomarkers]
        cluster_exp = pd.merge(all_tumor_cell_biomarkers, self.tumor_cell_cluster, on=["tumor", "cell"]) \
            .reset_index().set_index(["tumor", "cell", "cluster"])
        self.cluster_exp_quantile = cluster_exp.groupby(["tumor", "cluster"])[biomarkers].quantile(quantile_threshold)
        return self.cluster_exp_quantile

    def get_clusters_with_biomarker_expression(self, biomarkers):
        if type(biomarkers) is not list:
            biomarkers = list(biomarkers)
        tumor_cluster_ids = self.cluster_exp_quantile[(self.cluster_exp_quantile[biomarkers] > 0.0).any(axis=1)].index
        return tumor_cluster_ids







