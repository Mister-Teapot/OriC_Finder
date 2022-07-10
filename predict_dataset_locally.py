# Standard imports
import multiprocessing as mp
from functools import partial
import os, sys, csv, joblib
import numpy as np

# Set these before running
DATASET   = ['NC_002947', 'NC_002971']# Accessions that have been experimentally verified.
#     'NC_000964.1', 'NC_002947.1', 'NC_003272.1', 'NC_003869.1', 'NC_005090.1', 'NC_006461.1', 'NC_007633.1', 'NC_000913.1', 'NC_003047.1',
#     'NC_007604.1', 'NC_000962.1', 'NC_002696.1', 'NC_002971.1', 'NC_005363.1', 'NC_008255.1', 'NC_009850.1', 'NC_010546.1', 'NC_011916.1'
# ]
CSV_OUT_FOLDER = 'exp_csvs'
PARALLEL  = False
CPUS      = os.cpu_count() - 1
MAX_ORICS = 10 # Assumption: No more than 10 oriC for a single organism are predicted

EMAIL   = 'zoyavanmeel@gmail.com'
API_KEY = '795d705fb638507c9b2295c89cc64ee88108'
MODEL   = joblib.load('exp_train_model.pkl')

# Cluster path
sys.path.append('../OriC_Finder/')

from oriC_Finder import find_oriCs

def prep_prediction(accession, email, api_key, model, csv_path, max_oriCs):

    # Quick check to see if th FASTA has already been processed
    if os.path.exists(csv_path + '/' + accession + '.csv'):
        return

    # preferred_properties, all_oriCs = find_oriCs(sample_path)
    preferred_properties = find_oriCs(accession=accession, email=email, api_key=api_key, model=model)

    row = []
    RefSeq = preferred_properties['name'].split('.')[0] # RefSeq = RefSeq Accession Number, removing Version Number
    row.append(RefSeq) 
    row.append(preferred_properties['seq_size'])

    # Quality of Prediction processing
    row.append(preferred_properties['gc_conc'])

    # OriC processing
    for i in range(max_oriCs):
        row.append(preferred_properties['oriC_middles'][i]) if i < len(preferred_properties['oriC_middles']) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(preferred_properties['occurances'][i]) if i < len(preferred_properties['occurances']) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(preferred_properties['Z_occurances'][i]) if i < len(preferred_properties['Z_occurances']) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(preferred_properties['G_occurances'][i]) if i < len(preferred_properties['G_occurances']) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(preferred_properties['D_occurances'][i]) if i < len(preferred_properties['D_occurances']) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(preferred_properties['Prediction'][i]) if i < len(preferred_properties['Prediction']) else row.append(np.nan)

    with open(csv_path + '/' + RefSeq + '.csv', 'w') as fh:
        writer = csv.writer(fh)
        writer.writerow(row)
        fh.close()

if __name__ == '__main__':

    if not PARALLEL:
        for accession in DATASET:
            prep_prediction(accession, EMAIL, API_KEY, MODEL, CSV_OUT_FOLDER, MAX_ORICS)
    else:
        prepped_prediction = partial(prep_prediction, email=EMAIL, api_key=API_KEY, csv_path=CSV_OUT_FOLDER, max_oriCs=MAX_ORICS, model=MODEL)
        with mp.Pool(CPUS) as pool:
            pool.map(prepped_prediction, DATASET)