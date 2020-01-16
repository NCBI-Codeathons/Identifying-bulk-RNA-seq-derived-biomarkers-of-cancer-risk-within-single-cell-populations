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
        self.tumor_id = {tumor: tumor_id + 1 for tumor_id, tumor in enumerate(self.TUMORS)}
        self.id_tumor = {tumor_id + 1: tumor for tumor_id, tumor in enumerate(self.TUMORS)}


    def preprocess_data(self):
        self.tumor_ans = {}
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
            sc.pp.normalize_total(ad, target_sum=1e6)
            sc.pp.log1p(ad)
            ad.X = pd.DataFrame(ad.X, index=cells, columns=self.GENE_NAMES)
            self.tumor_ans[tumor] = ad

        for tumor in self.TUMORS:
            self.tumor_ans[tumor].X.index = pd.MultiIndex.from_tuples(self.tumor_ans[tumor].X.index.str.split("_", expand=True),
                                                                 names=["tumor", "cell"])
        self.all_tumor_df = pd.concat([self.tumor_ans[tumor].X for tumor in self.TUMORS])

    def counts_matrix_for_tumor(self, tumor):
        indices = self.TUMOR_TO_INDICES[tumor]
        with h5py.File(self.DATA_F, 'r') as f:
            counts = f['count'][indices]
        cells = list(np.array(self.CELLS)[indices])
        return counts, cells

    def get_unique_tumors(self):
        return np.unique(self.TUMORS)

    def add_clustering_results(self):
        pass





