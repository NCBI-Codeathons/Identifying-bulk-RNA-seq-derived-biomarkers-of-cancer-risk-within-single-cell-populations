from optparse import OptionParser
import pandas as pd
import json
import numpy as np
import h5py 

ROOT = '/Users/matthewbernstein/Development/single-cell-hackathon/data'
OUT_F = '/Users/matthewbernstein/Development/single-cell-hackathon/data/GSE76381.h5'
THE_FILE = '{}/GSE76381_EmbryoMoleculeCounts.cef.txt'.format(ROOT)

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
    print('Loading {}...'.format(THE_FILE))
    df = pd.read_csv(
        THE_FILE, 
        sep='\t', 
        header=None,
        index_col=0
    )
    cells = [
        x.encode('utf-8')
        for x in df.iloc[1,1:]
    ]
    cell_types = [
        x.encode('utf-8')
        for x in df.iloc[2,1:]
    ]
    df = df.iloc[5:,1:]
    df.columns = cells
    print(df)

    # Extract the genes
    gene_names = [
        x.encode('utf-8')
        for x in df.index
    ]

    assert df.shape[1] == len(cell_types)
    assert df.shape[1] == len(cells)
    assert df.shape[0] == len(gene_names)

    # Cast to a numpy array
    counts = np.array(df, dtype=np.int32).T

    print('Writing to H5 file...')
    with h5py.File(OUT_F, 'w') as f:
        f.create_dataset('count', data=counts)
        f.create_dataset('cell', data=np.array(cells))
        f.create_dataset('cell_type', data=np.array(cell_types))
        f.create_dataset('gene_name', data=np.array(gene_names))
    print('done.')


if __name__ == '__main__':
    main()
