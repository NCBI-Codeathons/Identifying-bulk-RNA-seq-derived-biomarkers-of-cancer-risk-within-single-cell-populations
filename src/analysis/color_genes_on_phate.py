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

sys.path.append('..')

from ..common import load_GSE103224 as load

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

def color_genes_on_phate(input_genes):
    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    for tumor in TUMORS:
        counts, cells = load.counts_matrix_for_tumor(tumor)
        ad = AnnData(
            X=counts, 
            obs=pd.DataFrame(data=cells, columns=['cell']),
            var=pd.DataFrame(
                index=load.GENE_NAMES, 
                data=load.GENE_NAMES, 
                columns=['gene_name']
            )
        )
        #print(ad_PJ048.var)

        sc.pp.normalize_total(ad, target_sum=1e6)
        sc.pp.log1p(ad)
       
        phate_operator = phate.PHATE(n_jobs=-2, random_state=1)
        X_phate = phate_operator.fit_transform(ad.X)
        ad.obs = pd.DataFrame(
            data=[[x,y,cell] for (x,y),cell in zip(X_phate, cells)], columns=['PHATE 1', 'PHATE 2', 'cell'])

        out_dir = join(tumor, 'cell_type_markers')
        os.system('mkdir -p {}'.format(out_dir))
        genes = input_genes
        _create_plots(ad, out_dir, genes)

        #out_dir = join(tumor, '6_survival_markers') 
        #os.system('mkdir -p {}'.format(out_dir))
        #genes = ['SLC16A3', 'MAP2K3', 'CD79B', 'IMPDH1', 'MPZL3', 'APOBR']
        #_create_plots(ad, out_dir, genes)

def _create_plots(ad, out_dir, genes):
    for gene in genes:
        # Cell type markers
        fig, ax = plt.subplots(1,1,figsize=(6,6))
        ax = sc.pl.scatter(ad, x='PHATE 1', y='PHATE 2', color=gene, ax=ax, legend_loc='none', show=False)
        ax.set_xticks([])
        ax.set_yticks([])
        l, b, w, h = fig.axes[-1].get_position().bounds
        ll, bb, ww, hh = fig.axes[0].get_position().bounds
        fig.axes[0].set_position([ll-0.05, bb, ww, hh])
        fig.axes[-1].set_position([ll+ww-0.05, b, w, h])
        #plt.tight_layout()
        fig.savefig(
            join(out_dir, '{}.png'.format(gene)),
            format='png',
            dpi=150
            #bbox_inches='tight'
        )

