# Standard imports
import multiprocessing as mp
from functools import partial
import csv
import os, sys
import numpy as np

# Set these before running
DATASET_NAME = 'refseq_23'
ON_CLUSTER   = False
PARALLEL     = True
MAX_ORICS    = 10 # Assumption: No more than 10 oriC for a single organism are predicted

# Self-made module
if ON_CLUSTER:
    # Cluster path
    sys.path.append('/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/OriC_Finder/')
    run_type = 'cluster'
else:
    # Local path
    sys.path.append('../OriC_Finder')
    run_type = 'local'
from oriC_Finder import find_oriCs, read_FASTA


def init_settings(run_type, dataset_name, parallel):
    '''Gets variables for processing based on whether the script is executed locally or on the cluster'''
    type_list = ['cluster', 'local']
    data_path, csv_path, cpus = None, None, None

    if run_type not in type_list:
        raise KeyError('Not a valid location')

    if run_type == 'cluster':
        data_path = '/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/' + dataset_name + '/bacteria'
        csv_path  = '/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/' + dataset_name + '/csvs'
        cpus      = 1

    if run_type == 'local':
        data_path = 'NCBI data prep/' + dataset_name + '/chromosomes_only'
        csv_path  = 'NCBI data prep/' + dataset_name + '/csvs'
        cpus      = os.cpu_count() - 1 if parallel else 1

    return data_path, csv_path, cpus


def prep_prediction(sample_path, csv_path, max_oriCs):
    # Quick check to see if th FASTA has already been processed
    with open(sample_path, 'r') as fh:
        accession = fh.readline().split(' ')[0][1:]
    if os.path.exists(csv_path + '/' + accession + '.csv'):
        return

    # preferred_properties, all_oriCs = find_oriCs(sample_path)
    preferred_properties = find_oriCs(sample_path)

    # Name processing
    name_list = preferred_properties['name'].split(' ')
    RefSeq = name_list[0]
    Organism = ' '.join(name_list[1:])
    row = []
    row.append(RefSeq) # RefSeq = RefSeq Accession Number
    row.append(Organism)
    row.append(preferred_properties['seq_size'])

    # Quality of Prediction processing
    row.append(preferred_properties['false_order']) # v1 & v2 &v3
    row.append(preferred_properties['gc_conc'])     # v3
    # row.append(preferred_properties['n_penalty']) # v1 & v2
    # row.append(preferred_properties['o_penalty']) # v1 & v2
    # row.append(preferred_properties['d_penalty']) # v1 & v2

    # OriC processing
    ## v1 & v2 & v3
    for i in range(max_oriCs):
        row.append(preferred_properties['oriC_middles'][i]) if i < len(preferred_properties['oriC_middles']) else row.append(np.nan)

    ## v3
    for i in range(max_oriCs):
        row.append(preferred_properties['occurances'][i]) if i < len(preferred_properties['occurances']) else row.append(np.nan)
    
    ## v1 & v2
    # for i in range(max_oriCs):
    #     row.append(preferred_properties['oriC_edges'][i]) if i < len(preferred_properties['oriC_edges']) else row.append(np.nan)
    # for i in range(max_oriCs):
    #     row.append(all_oriCs['xy_oriCs'][i]) if i < len(all_oriCs['xy_oriCs']) else row.append(np.nan)
    # for i in range(max_oriCs):
    #     row.append(all_oriCs['xgc_oriCs'][i]) if i < len(all_oriCs['xgc_oriCs']) else row.append(np.nan)
    # for i in range(max_oriCs):
    #     row.append(all_oriCs['ygc_oriCs'][i]) if i < len(all_oriCs['ygc_oriCs']) else row.append(np.nan)

    with open(csv_path + '/' + RefSeq + '.csv', 'w') as fh:
        writer = csv.writer(fh)
        writer.writerow(row)
        fh.close()

if __name__ == '__main__':
    data_path, csv_path, cpus = init_settings(run_type, DATASET_NAME, PARALLEL)
    samples = os.listdir( data_path )
    sample_paths = [data_path + '/' + fasta for fasta in samples]
    del samples
    if not PARALLEL:
        for sample_path in sample_paths:
            prep_prediction(sample_path, csv_path, MAX_ORICS)
    else:    
        with mp.Pool(cpus) as pool:
            prepped_prediction = partial(prep_prediction, csv_path=csv_path, max_oriCs=MAX_ORICS)
            pool.map(prepped_prediction, sample_paths)
