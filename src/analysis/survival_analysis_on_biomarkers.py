import matplotlib as mpl
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

import sklearn
from sklearn.neighbors import BallTree
import scipy
from scipy.spatial import distance
from scipy.stats import zscore
import sys
import numpy as np
import pandas as pd
import lifelines
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

sys.path.append('/Users/matthewbernstein/Development/single-cell-hackathon/Identifying-bulk-RNA-seq-derived-biomarkers-of-cancer-risk-within-single-cell-populations/src/common')

import load_TCGA_GBM
import load_GSE103224
import normalize
import slice_matrix

BIOMARKERS_F = '/Users/matthewbernstein/Development/single-cell-hackathon/Identifying-bulk-RNA-seq-derived-biomarkers-of-cancer-risk-within-single-cell-populations/data/MK_genes_TUMORS_integrated.csv'

def main():
    df = pd.read_csv(BIOMARKERS_F)
    MARKERS = list(df['gene'])

    # Retrieve counts matrix
    X, samples = load_TCGA_GBM.counts_matrix()

    # Normalize
    X = normalize.log_cpm(X)

    X = slice_matrix.keep_genes(
        X, 
        load_TCGA_GBM.GENE_NAMES, 
        MARKERS
    )

    # Compute x-scores
    X = zscore(X)

    # Make sure the metadata and matrix are aligned
    assert tuple(load_TCGA_GBM.SAMPLE_IDS) == tuple(load_TCGA_GBM.TCGA_GBM_META['sample_id']) 

    if X.shape[1] > 1:
        scores = np.sum(X, axis=1)
    else:
        scores = np.squeeze(X, axis=1)
        
    up_q = 0.6 # Top quantile
    low_q = 0.4 # Bottom quantile
    high = np.quantile(scores, up_q)
    low = np.quantile(scores, low_q)

    df = load_TCGA_GBM.TCGA_GBM_META
    df_scores = pd.DataFrame(
        data=[
            (sample, score) 
            for sample,score in zip(load_TCGA_GBM.SAMPLE_IDS, scores)
        ],
        columns=['sample_id', 'score']
    )
    df_scores = df_scores.set_index('sample_id')
    df = df.set_index('sample_id')
    df = df.join(df_scores, on='sample_id', how='left')


    df_high = df.loc[df['score'] > high]
    df_low = df.loc[df['score'] <= low]


    r = logrank_test(
        df_high['time'],
        df_low['time'],
        event_observed_A=df_high['censor'],
        event_observed_B=df_low['censor']
    )
    logrank_pval = r.p_value
    print('log-rank p-value: {}'.format(logrank_pval))

    kmf = KaplanMeierFitter()
    kmf.fit(
        df_high['time'],
        df_high['censor'],
        label='>{}th Quantile (n={})'.format(int(up_q * 100), len(df_high))
    )
    ax = kmf.plot(ci_show=True, show_censors=True, color='hotpink', ls='-')
    kmf.fit(
        df_low['time'],
        df_low['censor'],
        label='<{}th Quantile (n={})'.format(int(low_q * 100), len(df_low))
    )
    ax = kmf.plot(ci_show=True, show_censors=True, ax=ax, color='mediumblue', ls='-')
    ax.set_xlim((0,730))
    plt.show()

if __name__ == '__main__':
    main()
