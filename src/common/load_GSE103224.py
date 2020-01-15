import h5py
import numpy as np
from collections import defaultdict

DATA_F = '/Users/matthewbernstein/Development/single-cell-hackathon/data/GSE103224.h5'

# Load cell names and tumor names
with h5py.File(DATA_F, 'r') as f:
    CELLS = [
        str(x)[2:-1]
        for x in f['cell'][:]
    ]
    TUMORS = [
        str(x)[2:-1]
        for x in f['tumor'][:]
    ]
    GENE_IDS = [
        str(x)[2:-1]
        for x in f['gene_id'][:]
    ]
    GENE_NAMES = [
        str(x)[2:-1]
        for x in f['gene_name'][:]
    ]

# Map each cell to its index in the data matrix
CELL_TO_INDEX = {
    cell: index
    for index, cell in enumerate(CELLS)
}

# Map each tumor to its indices in the data matrix
TUMOR_TO_INDICES = defaultdict(lambda: [])
for index, tumor in enumerate(TUMORS):
    TUMOR_TO_INDICES[tumor].append(index)
TUMOR_TO_INDICES = dict(TUMOR_TO_INDICES)


def counts_matrix_for_tumor(tumor):
    """
    Retrieve the counts matrix for a given tumor.

    Args:
        tumor: the tumor ID
    Returns:
        counts: matrix of counts for tumor
        cells: cell names correponding to the
            rows of counts
    """
    indices = TUMOR_TO_INDICES[tumor]
    with h5py.File(DATA_F, 'r') as f:
        counts = f['count'][indices]
    cells = list(np.array(CELLS)[indices])
    return counts, cells


def counts_matrix_for_cells(cells):
    """
    Retrieve the counts matrix for a specific
    set of cells.

    Args:
        cells: list of N cells
    Returns:
        NxM matrix of counts for the N
        cells and M genes. Rows are returned
        in same order as input list of cells.
    """
    indices = [
        CELL_TO_INDEX[cell]
        for cell in cells
    ]
    with h5py.File(DATA_F, 'r') as f:
        counts = f['count'][indices]
    return counts


def main():
    cells = ['PJ016_4', 'PJ016_5', 'PJ016_6']
    #print(np.sum(counts_matrix_for_cells(cells), axis=1))

if __name__ == '__main__':
    main()
