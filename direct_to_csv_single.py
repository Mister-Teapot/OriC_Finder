# Standard imports
import multiprocessing as mp
import csv
import os, sys
import numpy as np

# Set these before running
DATASET_NAME = 'refseq_3k_set'
ON_CLUSTER   = True
PARALLEL     = False

# Self-made module
if ON_CLUSTER:
    # Cluster path
    sys.path.append('/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/OriC_Finder/')
    run_type = 'cluster'
else:
    # Local path
    sys.path.append('../OriC_Finder')
    run_type = 'local'
from oriC_Finder import find_oriCs


# df['RefSeq'] = df['RefSeq'].str.extract(r'([^.]*)')

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
        print('mp.cpu_count() found :', mp.cpu_count())

    if run_type == 'local':
        path   = os.path.join('NCBI data prep', dataset_name, 'chromosomes_only')
        to_csv = os.path.join('NCBI data prep', dataset_name)
        cpus   = mp.cpu_count() - 1 if parallel else 1

    return path, to_csv, cpus

if __name__ == '__main__':
    path, to_csv, _ = get_standard_vars(run_type, DATASET_NAME, PARALLEL)
    samples = os.listdir( path )
    sample_paths = [os.path.join(path, fasta) for fasta in samples]

    fieldnames = [
        'RefSeq',
        'Organism',
        'Sequence_length',
        'Penalties',
        'False_order'
    ]

    max_orics = 10 # Assumption: No more than 10 oriC for a single organism are predicted
    for i in range(max_orics):
        fieldnames.append(f'oriC_edges_{i}' )
        fieldnames.append(f'oriC_middles_{i}')

    with open('get_whatever_i_can.csv', 'w') as fh:
        writer_object = csv.writer(fh)
        writer_object.writerow(fieldnames)
        fh.close()

    for sample in sample_paths:
        properties = find_oriCs(sample)
        # Name processing
        name_list = properties['name'].split(' ')
        RefSeq = name_list[0]
        Organism = ' '.join(name_list[1:])
        row = []
        row.append(RefSeq) # RefSeq = RefSeq Accession Number
        row.append(Organism)
        row.append(properties['seq_size'])

        # Quality of Prediction processing
        row.append(properties['nod_penalties'])
        row.append(properties['false_order'])

        # OriC processing
        for i in range(max_orics):
            if i < len(properties['oriC_middles']):
                row.append(properties['oriC_edges'][i])
                row.append(properties['oriC_middles'][i])
            else:
                row.append(np.nan)
                row.append(np.nan)

        with open('get_whatever_i_can.csv', 'a', newline='') as fh:
            writer_object = csv.writer(fh)
            writer_object.writerow(row)
            fh.close()
