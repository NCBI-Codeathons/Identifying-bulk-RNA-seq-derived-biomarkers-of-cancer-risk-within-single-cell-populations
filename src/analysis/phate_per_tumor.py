import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import sys
import os 
from os.path import join
import subprocess
from anndata import AnnData
import phate
from optparse import OptionParser

# TODO THIS NEEDS TO BE CHANGED
sys.path.append('src/common')

DATA_F = 'GSE103224.h5'

from load_h5 import H5COUNTS

TUMORS = [
    'PJ016',
    'PJ018',
    'PJ048',
    'PJ030',
    'PJ025',
    'PJ035',
    'PJ017',
    'PJ032'
]

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    (options, args) = parser.parse_args()
    
    data_root = args[0]
    cancer_biomarker = args[1]
    cell_type_biomarkers = args[2].split(',')
    out_dir = options.out_dir

    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    data = H5COUNTS(join(data_root, DATA_F))
    data.preprocess_data()
    data.add_clustering_results(path=join(data_root, 'interim/'))

    for tumor in TUMORS:
        ad = data.tumor_to_ad[tumor]

        obs_filt = ad.obs.loc[
            ad.obs['cluster'].notnull()
        ]
        indices = [int(x) for x in obs_filt.index]
        X_filt = ad.X.iloc[indices]
        X_filt = X_filt.set_index(obs_filt.index)
        ad_filt = AnnData(
            X=X_filt, 
            obs=obs_filt,
            var=ad.var
        )

        phate_operator = phate.PHATE(n_jobs=-2, random_state=1)
        X_phate = phate_operator.fit_transform(ad_filt.X)
        ad_filt.obs = pd.DataFrame(
            data=[
                [x, y, cluster] 
                for (x,y), cluster in 
                zip(X_phate, ad_filt.obs['cluster'])
            ], 
            columns=[
                'PHATE 1', 
                'PHATE 2', 
                'cluster'
            ]
        )

        # Color points by cluster
        fig, ax = plt.subplots(1,1,figsize=(8,6))
        ax = sc.pl.scatter(
            ad_filt, 
            x='PHATE 1', 
            y='PHATE 2', 
            color='cluster', 
            ax=ax, 
            legend_loc='right margin', 
            show=False
        )
        ax.set_xticks([])
        ax.set_yticks([])
        l, b, w, h = fig.axes[-1].get_position().bounds
        ll, bb, ww, hh = fig.axes[0].get_position().bounds
        plt.tight_layout()
        fig.savefig(
            join(out_dir, '{}_color_by_cluster.png'.format(tumor)),
            format='png',
            dpi=150
            #bbox_inches='tight'
        )

        # Color by genes
        genes = [cancer_biomarker] + cell_type_biomarkers
        for gene in genes:
            fig, ax = plt.subplots(1,1,figsize=(8,6))
            ax = sc.pl.scatter(
                ad_filt,
                x='PHATE 1',
                y='PHATE 2',
                color=gene,
                ax=ax,
                legend_loc='right margin',
                show=False
            )
            ax.set_xticks([])
            ax.set_yticks([])
            l, b, w, h = fig.axes[-1].get_position().bounds
            ll, bb, ww, hh = fig.axes[0].get_position().bounds
            plt.tight_layout()
            fig.savefig(
                join(out_dir, '{}_color_by_{}.png'.format(tumor, gene)),
                format='png',
                dpi=150
                #bbox_inches='tight'
            )

    
if __name__ == '__main__':
    main()

