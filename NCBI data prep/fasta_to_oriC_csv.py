# Standard imports
import os, sys
import pandas as pd
import numpy as np

# Local imports
sys.path.append('../OriC_Finder')
from oriC_Finder import find_oriCs

# Pandas printing options
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 30)


def database_oriC_prediction(path, to_csv=None):
    """
    Predict oriCs for database at given path.
    Compiles everything in a CSV-file in the given folder location.
    Returns pandas.DataFrame
    """

    sample_list = os.listdir( path )
    df_dict = {
        'RefSeq'         : [],
        'Organism'       : [],
        'Sequence_length': [],
        'Penalties'      : [],
        'False_order'    : []
    }
    max_orics = 10 # Assumption: No more than 10 oriC for a single organism are predicted
    for i in range(max_orics):
        df_dict.update( {f'oriC_edges_{i}': []} )
        df_dict.update( {f'oriC_middles_{i}': []} )

    for fasta in sample_list:
        properties = find_oriCs( os.path.join(path, fasta) )
        
        # Name processing
        name_list = properties['name'].split(' ')
        RefSeq = name_list[0]
        Organism = ' '.join(name_list[1:])

        df_dict['RefSeq'].append(RefSeq) # RefSeq = RefSeq Accession Number
        df_dict['Organism'].append(Organism)
        df_dict['Sequence_length'].append(properties['seq_size'])

        # Quality of Prediction processing
        df_dict['Penalties'].append(properties['nod_penalties'])
        df_dict['False_order'].append(properties['false_order'])

        # OriC processing
        for i in range(max_orics):
            if i < len(properties['oriC_middles']):
                df_dict[f'oriC_edges_{i}'].append(properties['oriC_edges'][i])
                df_dict[f'oriC_middles_{i}'].append(properties['oriC_middles'][i])
            else:
                df_dict[f'oriC_edges_{i}'].append(np.nan)
                df_dict[f'oriC_middles_{i}'].append(np.nan)

    df = pd.DataFrame.from_dict(df_dict)
    df.dropna(how='all', axis=1, inplace=True)

    # Removing version numbers for now
    df['RefSeq'] = df['RefSeq'].str.extract(r'([^.]*)')

    if to_csv is not None:
        df.to_csv(to_csv + f'/NCBI_oriC_{df.shape[0]}_improved.csv', index=False)
    return df

if __name__ == '__main__':
    _ = database_oriC_prediction('NCBI data prep/refseq_15/chromosomes_only', to_csv='NCBI data prep/refseq_15')