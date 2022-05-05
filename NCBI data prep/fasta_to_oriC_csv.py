# Standard imports
import multiprocessing as mp
import os, sys
import pandas as pd
import numpy as np

# Set these before running
DATASET_NAME = 'refseq_15'
ON_CLUSTER   = False
PARALLEL     = True

# Self-made module
if ON_CLUSTER:
    # Cluster path
    sys.path.append('/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/OriC_Finder/')
else:
    # Local path
    sys.path.append('../OriC_Finder')
from oriC_Finder import find_oriCs

# Pandas printing options
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 30)


def database_oriC_prediction(properties_list, to_csv=None):
    """
    Predict oriCs for database at given path.
    Compiles everything in a CSV-file in the given folder location.
    Returns pandas.DataFrame
    """

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

    for properties in properties_list:    
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
        df.to_csv(to_csv + f'/NCBI_oriC_{df.shape[0]}.csv', index=False)
    return df


def get_standard_vars(run_type, dataset_name, parallel):
    '''Gets variables for processing based on whether the script is executed locally or on the cluster'''
    type_list = ['cluster', 'local']
    path      = None
    to_csv    = None
    cpus      = None

    if run_type not in type_list:
        raise KeyError('Not a valid location')

    if run_type == 'cluster':
        path   = '/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/' + dataset_name + '/bacteria'
        to_csv = '/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/' + dataset_name
        cpus   = mp.cpu_count() if parallel else 1

    if run_type == 'local':
        path   = os.path.join('NCBI data prep', dataset_name, 'chromosomes_only')
        to_csv = os.path.join('NCBI data prep', dataset_name)
        cpus   = mp.cpu_count() - 1 if parallel else 1

    return path, to_csv, cpus

if __name__ == '__main__':
    path, to_csv, cpus = get_standard_vars('cluster', DATASET_NAME, PARALLEL)
    samples = os.listdir( path )
    sample_paths = [os.path.join(path, fasta) for fasta in samples]

    with mp.Pool(cpus) as p:
        prop_list = p.map(find_oriCs, sample_paths )
    
    database_oriC_prediction(properties_list=prop_list, to_csv=to_csv)
    