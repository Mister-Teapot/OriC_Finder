# Standard imports
import multiprocessing as mp
from functools import partial
import csv
import os, sys
import numpy as np

# Set these before running
DATASET_NAME = 'not_predicted.txt'
PARALLEL     = True
CPUS         = 3
MAX_ORICS    = 10 # Assumption: No more than 10 oriC for a single organism are predicted

EMAIL   = 'zoyavanmeel@gmail.com'
API_KEY = '795d705fb638507c9b2295c89cc64ee88108'

# Cluster path
sys.path.append('../OriC_Finder/')

from oriC_Finder import find_oriCs

def prep_prediction(accession, email, api_key, csv_path, max_oriCs):

    # Quick check to see if th FASTA has already been processed
    if os.path.exists(csv_path + '/' + accession + '.csv'):
        return

    # preferred_properties, all_oriCs = find_oriCs(sample_path)
    preferred_properties = find_oriCs(accession=accession, email=email, api_key=api_key)

    row = []
    RefSeq = preferred_properties['name'].split('.')[0] # RefSeq = RefSeq Accession Number, removing Version Number
    row.append(RefSeq) 
    row.append(preferred_properties['seq_size'])

    # Quality of Prediction processing
    row.append(preferred_properties['false_order'])
    row.append(preferred_properties['gc_conc'])

    # OriC processing
    for i in range(max_oriCs):
        row.append(preferred_properties['oriC_middles'][i]) if i < len(preferred_properties['oriC_middles']) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(preferred_properties['occurances'][i]) if i < len(preferred_properties['occurances']) else row.append(np.nan)

    with open(csv_path + '/' + RefSeq + '.csv', 'w') as fh:
        writer = csv.writer(fh)
        writer.writerow(row)
        fh.close()

if __name__ == '__main__':
    csv_path       = DATASET_NAME.strip('.txt') + '_csvs'

    with open(DATASET_NAME) as fh:
        not_predicted = fh.read().split('\n')

    if not PARALLEL:
        for accession in not_predicted:
            prep_prediction(accession, EMAIL, API_KEY, csv_path, MAX_ORICS)
    else:
        prepped_prediction = partial(prep_prediction, email=EMAIL, api_key=API_KEY, csv_path=csv_path, max_oriCs=MAX_ORICS)
        with mp.Pool(CPUS) as pool:
            pool.map(prepped_prediction, not_predicted)