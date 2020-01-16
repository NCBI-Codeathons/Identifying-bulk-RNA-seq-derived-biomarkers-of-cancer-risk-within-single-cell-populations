import h5py
import numpy as np
from collections import defaultdict
import pandas as pd

EXPRESSION_F = '/Users/matthewbernstein/Development/single-cell-hackathon/data/TCGA_GBM.h5'
METADATA_F = '/Users/matthewbernstein/Development/single-cell-hackathon/data/TCGA_GBM_metadata.tsv'

TCGA_GBM_META = pd.read_csv(METADATA_F, sep='\t')

# Load sample ID's and gene ID's
with h5py.File(EXPRESSION_F, 'r') as f:
    GENE_IDS = [
        str(x)[2:-1]
        for x in f['gene_id'][:]
    ]
    GENE_NAMES = [
        str(x)[2:-1]
        for x in f['gene_name'][:]
    ]
    SAMPLE_IDS = [
        str(x)[2:-1]
        for x in f['sample_id'][:]
    ]

# Map each sample to its index in the data matrix
SAMPLE_TO_INDEX = {
    sample: index
    for index, sample in enumerate(SAMPLE_IDS)
}

def counts_matrix():
    with h5py.File(EXPRESSION_F, 'r') as f:
        counts = f['count'][:]
    return counts, SAMPLE_IDS

def counts_matrix_for_samples(samples=None):
    """
    Retrieve the counts matrix for a specific
    set of samples.

    Args:
        samples: list of N samples
    Returns:
        NxM matrix of counts for the N
        samples and M genes. Rows are returned
        in same order as input list of samples.
    """
    indices = [
        SAMPLE_TO_INDEX[sample]
        for sample in samples
    ]
    with h5py.File(EXPRESSION_F, 'r') as f:
        counts = f['count'][indices]
    return counts

def main():
    samples = [
        'fe7d8f81-f39f-4401-b3e9-cfde245d6b93',
        '53576d9b-d84c-40ac-8a82-ac6bc81180a2',
        '2e54b242-9fda-4cb6-aeb7-4c3a56f8871a'
    ]
    #print(counts_matrix_for_samples(samples))
    print(np.sum(counts_matrix_for_samples(samples), axis=1))

if __name__ == '__main__':
    main()
