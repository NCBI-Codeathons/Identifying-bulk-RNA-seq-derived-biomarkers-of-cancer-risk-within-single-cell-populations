######################################
#   Basic normalization functions
######################################
import numpy as np

def keep_genes(X, genes, keep_genes):
    gene_to_index = {
        gene: index
        for index, gene in enumerate(genes)
    }
    indices = [
        gene_to_index[gene]
        for gene in keep_genes
        if gene in gene_to_index
    ]
    return X[:,indices]
