from optparse import OptionParser
import pandas as pd
import json
import numpy as np
import h5py 

ROOT = '/Users/matthewbernstein/Development/single-cell-hackathon/data/GSE103224_RAW'
OUT_F = '/Users/matthewbernstein/Development/single-cell-hackathon/data/GSE103224.h5'
TUMOR_TO_FILE = {
    'PJ016': '{}/GSM2758471_PJ016.filtered.matrix.txt'.format(ROOT),
    'PJ018': '{}/GSM2758473_PJ018.filtered.matrix.txt'.format(ROOT),
    'PJ030': '{}/GSM2758475_PJ030.filtered.matrix.txt'.format(ROOT),
    'PJ035': '{}/GSM2758477_PJ035.filtered.matrix.txt'.format(ROOT),
    'PJ017': '{}/GSM2758472_PJ017.filtered.matrix.txt'.format(ROOT),
    'PJ025': '{}/GSM2758474_PJ025.filtered.matrix.txt'.format(ROOT),
    'PJ032': '{}/GSM2758476_PJ032.filtered.matrix.txt'.format(ROOT),
    'PJ048': '{}/GSM2940098_PJ048.filtered.matrix.txt'.format(ROOT)
}

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    #parser.add_option("-a", "--a_descrip", action="store_true", help="This is a flat")
    #parser.add_option("-b", "--b_descrip", help="This is an argument")
    (options, args) = parser.parse_args()

    print('Finding genes common to all datasets...')
    the_genes = None
    counts = None
    cells = []
    tumors = []
    for tumor, tumor_f in TUMOR_TO_FILE.items():
        print('Loading {}...'.format(tumor_f))
        df = pd.read_csv(
            tumor_f, 
            sep='\t', 
            header=None, 
            index_col=0
        )

        # Extract the genes
        genes = list(df.index)
        if the_genes is None:
            the_genes = list(genes)
            gene_names = list(df.iloc[:,0])
        assert frozenset(genes) == frozenset(the_genes)

        # Drop the gene-names column
        keep_cols = df.columns[1:]
        df = df[keep_cols]

        # Re-order rows according to the global gene 
        # ordering
        df = df.loc[the_genes]

        # Cast to a numpy array
        curr_counts = np.array(df).T
        if counts is None:
            counts = curr_counts
        else:
            counts = np.concatenate([counts, curr_counts])
        print('Current shape of the expression matrix: {}'.format(counts.shape))

        cells += [
            '{}_{}'.format(tumor, i).encode('utf-8')
            for i in np.arange(curr_counts.shape[0])
        ]
        tumors += [
            tumor.encode('utf-8')
            for i in np.arange(curr_counts.shape[0])
        ]


    the_genes = [
        x.encode('utf-8')
        for x in the_genes
    ]
    gene_names = [
        x.encode('utf-8')
        for x in gene_names
    ]
    print('done.')

    print('Writing to H5 file...')
    with h5py.File(OUT_F, 'w') as f:
        f.create_dataset('count', data=counts)
        f.create_dataset('cell', data=np.array(cells))
        f.create_dataset('tumor', data=np.array(tumors))
        f.create_dataset('gene_id', data=np.array(the_genes))
        f.create_dataset('gene_name', data=np.array(gene_names))
    print('done.')


        #print(genes[:50])
        #assert frozenset([x.split('.')[0] for x in genes]) == frozenset([x.split('.')[0] for x in the_genes])
        #print(len(genes))
       # print(genes[:50])

    #genes = df['Unnamed: 0']
    #print(genes)

    #df = df[keep_cols]
    #df.rename(
    #    index={
    #        i: gene
    #        for i,gene in enumerate(genes)
    #    },
    #    inplace=True
    #)
    #df.to_csv('../GSE57872_processed.tsv', sep='\t')

if __name__ == '__main__':
    main()
