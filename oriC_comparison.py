import os
# Set cwd to location of this script
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

from Bio import pairwise2 # Apperently it's better to use Align, but I can't get it to work for now
import pandas as pd
import re

Z_curve_df = pd.read_csv( os.path.join('NCBI_data_prep', 'NCBI_oriC_20_improved.csv') )
DoriC_df   = pd.read_csv( os.path.join('DoriC_data_prep', 'DoriC_oriC_split.csv') )

checked = 0
partial_oriC = []
for DoriC_idx, accession in enumerate(DoriC_df['RefSeq']):
    for Z_curve_idx, accession in enumerate(list(Z_curve_df['RefSeq'])):
        all_idx = [x for x in range(Z_curve_df['oriC_i'][Z_curve_idx], Z_curve_df['oriC_f'][Z_curve_idx]+1)]

        if DoriC_df['oriC_i'][DoriC_idx] in all_idx or DoriC_df['oriC_f'][DoriC_idx] in all_idx:
            checked += 1
            partial_oriC.append(accession)

partial_oriC = set(partial_oriC)
print(f"{checked} combinations were in checked")
print(f"{len(partial_oriC)} had partial oriC overlap.\n")
for i in partial_oriC:
    print(i)