# Standard imports
import multiprocessing as mp
from functools import partial
import csv
import os, sys
import numpy as np

# Set these before running
DATASET_NAME = 'refseq_DoriC_accessions_set'
PARALLEL     = True
CPUS         = 10
MAX_ORICS    = 10 # Assumption: No more than 10 oriC for a single organism are predicted

# Cluster path
sys.path.append('/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/OriC_Finder/')

from oriC_Finder import find_oriCs

def prep_prediction(sample_path, csv_path, max_oriCs):
    gene_info_path = sample_path[0]
    fasta_path = sample_path[1]

    # Quick check to see if th FASTA has already been processed
    with open(fasta_path, 'r') as fh:
        accession = fh.readline().split(' ')[0][1:]
    if os.path.exists(csv_path + '/' + accession + '.csv'):
        return

    # preferred_properties, all_oriCs = find_oriCs(sample_path)
    preferred_properties = find_oriCs(genome_fasta=fasta_path, genes_fasta=gene_info_path)

    row = []
    RefSeq = preferred_properties['name'][:-2] # RefSeq = RefSeq Accession Number, removing Version Number
    row.append(RefSeq) 
    row.append(preferred_properties['seq_size'])

    # Quality of Prediction processing
    row.append(preferred_properties['false_order'])
    row.append(preferred_properties['gc_conc'])

    # OriC processing
    for i in range(max_oriCs):
        row.append(preferred_properties['oriC_middles'][i]) if i < len(preferred_properties['oriC_middles']) else row.append(np.nan)
    ## v3
    for i in range(max_oriCs):
        row.append(preferred_properties['occurances'][i]) if i < len(preferred_properties['occurances']) else row.append(np.nan)

    with open(csv_path + '/' + RefSeq + '.csv', 'w') as fh:
        writer = csv.writer(fh)
        writer.writerow(row)
        fh.close()

if __name__ == '__main__':
    fastas_path    = '/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/' + DATASET_NAME + '/bacteria'
    gene_info_path = '/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/' + DATASET_NAME + '/gene_info_files'
    csv_path       = '/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/' + DATASET_NAME + '/csvs'

    samples = os.listdir( gene_info_path )
    sample_paths = [(gene_info_path + '/' + sample, fastas_path + '/' + sample) for sample in samples]
    del samples

    if not PARALLEL:
        for path in sample_paths:
            prep_prediction(path, csv_path, MAX_ORICS)
    else:
        prepped_prediction = partial(prep_prediction, csv_path=csv_path, max_oriCs=MAX_ORICS)
        with mp.Pool(CPUS) as pool:
            pool.map(prepped_prediction, sample_paths)