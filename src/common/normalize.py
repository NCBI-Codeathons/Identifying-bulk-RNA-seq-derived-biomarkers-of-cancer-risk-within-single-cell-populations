######################################
#   Basic normalization functions
######################################
import numpy as np

def cpm(X, const=10e6):
    X_norm = np.array([
        x/sum(x)
        for x in X
    ])
    X_norm *= const
    return X_norm


def log_cpm(X):
    X_cpm = cpm(X)
    return np.log(X_cpm+1)


def restrict_common_genes(X1, X2, genes1, genes2):
    """
    Given two expression matrices and two lists of gene
    identifiers corresponding to the two columns of the
    matrices, re-configure the matrices to keep only
    columns for genes that are common between the two
    datasets.
    """
    gene_to_index1 = {
        gene: index
        for index, gene in enumerate(genes1)
    }
    gene_to_index2 = {
        gene: index
        for index, gene in enumerate(genes2)
    }
    common_genes = sorted(
        set(genes1) & set(genes2)
    )
    indices1 = [
        gene_to_index1[gene]
        for gene in common_genes
    ]
    indices2 = [
        gene_to_index2[gene]
        for gene in common_genes
    ]
    X1_new = X1[:,indices1]
    X2_new = X2[:,indices2]
    return X1_new, X2_new, common_genes
