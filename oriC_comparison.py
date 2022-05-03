import os
from ast import literal_eval
import pandas as pd
import re

# Set cwd to location of this script
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

from Bio import pairwise2 # Apperently it's better to use Align, but I can't get it to work for now

# Pandas printing options
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 25)


# Read CSVs
Z_curve_df = pd.read_csv( 'NCBI data prep/NCBI_oriC_15_improved.csv' )
DoriC_df   = pd.read_csv( 'DoriC data prep/DoriC_oriC_seperate_entries.csv' )

# Setting columns to proper types
Z_curve_df = Z_curve_df.astype({'Sequence_length':'float', 'False_order':'bool'})
Z_curve_df['Penalties'] = Z_curve_df['Penalties'].map(literal_eval)

oriC_cols = [i for i in Z_curve_df.columns if 'oriC' in i]
most_oriC = int(sorted(oriC_cols, key=lambda x:x[-1], reverse=True)[0][-1]) + 1
for i in range(most_oriC):
    Z_curve_df = Z_curve_df.astype({f'oriC_middles_{i}':'float'})

    Z_curve_df[f'oriC_edges_{i}'].fillna('()', inplace=True)
    Z_curve_df[f'oriC_edges_{i}'] = Z_curve_df[f'oriC_edges_{i}'].map(literal_eval)

# All DoriC columns have been interpreted properly
# print(DoriC_df.dtypes)

# checked = 0
# partial_oriC = []
# for DoriC_idx, accession in enumerate(DoriC_df['RefSeq']):
#     for Z_curve_idx, accession in enumerate(list(Z_curve_df['RefSeq'])):
#         all_idx = [x for x in range(Z_curve_df['oriC_i'][Z_curve_idx], Z_curve_df['oriC_f'][Z_curve_idx]+1)]

#         if DoriC_df['oriC_i'][DoriC_idx] in all_idx or DoriC_df['oriC_f'][DoriC_idx] in all_idx:
#             checked += 1
#             partial_oriC.append(accession)

# partial_oriC = set(partial_oriC)
# print(f"{checked} combinations were in checked")
# print(f"{len(partial_oriC)} had partial oriC overlap.\n")
# for i in partial_oriC:
#     print(i)