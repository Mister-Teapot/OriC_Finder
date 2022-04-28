import re
from matplotlib.pyplot import axis
import numpy as np

import pandas as pd
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 15)

import os
# Set cwd to location of this script
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

'''
oriFinder will always produce a ori 250 nucleotides upstream of the skewplot
minimum and 700 nucleotides downstream. This may be possible to further optimise
using DnaA-box regions.The problem there is that the DnaA-box region of many
organisms are unknown.

This script aims to evaluate the performance of oriFinder against the oriC
regions predicted and compiled in DoriC.
'''


# Tianjin University BioInformatics Centre (tubic) a.k.a. the people who made DoriC
# Copied Katelyn's processing of the csv with a few modifications
raw_oriC = pd.read_csv('tubic_bacteria.csv')

# Remove duplicates
raw_oriC[raw_oriC['Refseq'].duplicated(keep = False) == True]
raw_oriC.rename(columns = {'Refseq' : 'RefSeq'}, inplace = True)

# Extract useful columns
raw_oriC = raw_oriC[['RefSeq', 'Organism', 'Lineage', 'Location of replication genes', 'Location of replication origin', 'OriC sequence']].copy()

# Give beginning and end point of each oriC and DnaA a separate column
raw_oriC['Location of replication origin'] = raw_oriC['Location of replication origin'].str.split(r'\.\.')
raw_oriC['oriC_i'] = raw_oriC['Location of replication origin'].str[0]
raw_oriC['oriC_f'] = raw_oriC['Location of replication origin'].str[1]
raw_oriC['oriC_i'] = raw_oriC['oriC_i'].str.extract(r'(\d+)')
raw_oriC['oriC_f'] = raw_oriC['oriC_f'].str.extract(r'(\d+)')
del raw_oriC['Location of replication origin']

# If an organisms has multiple DnaA boxes, they do not get separate entries, so don't forget to check for them!
split_DnaA = raw_oriC['Location of replication genes'].str.split(r',', expand=True)

raw_oriC = pd.concat( [raw_oriC, split_DnaA] , axis=1 )
raw_oriC.drop('Location of replication genes', axis=1, inplace=True)
for i in range(split_DnaA.shape[1]):
    # Split start and end of each DnaA box
    raw_oriC[i] = raw_oriC[i].str.split(r'\.\.')

    # Give start and end their own column
    raw_oriC[f'DnaA_{i}_i'] = raw_oriC[i].str[0]
    raw_oriC[f'DnaA_{i}_f'] = raw_oriC[i].str[1]

    # Remove anything that is not a number
    raw_oriC[f'DnaA_{i}_i'] = raw_oriC[f'DnaA_{i}_i'].str.extract(r'(\d+)')
    raw_oriC[f'DnaA_{i}_f'] = raw_oriC[f'DnaA_{i}_f'].str.extract(r'(\d+)')

    raw_oriC = raw_oriC.astype( {f'DnaA_{i}_i':'float', f'DnaA_{i}_f':'float'} )

    # Remove the unsplit column
    raw_oriC.drop(i, axis=1, inplace=True)

raw_oriC = raw_oriC.astype( {f'oriC_i':'float', 'oriC_f':'float'} )

# Remove RefSeq versions in accension number (might keep them later, since DoriC is from 2018)
raw_oriC['RefSeq'] = raw_oriC['RefSeq'].str.extract(r'([^.]*)')
raw_oriC['Lineage'] = raw_oriC['Lineage'].str.rstrip('.')

# NOT removing multichromosomal organisms -> Organisms with multiple chromosomes already have separate entries

# Renaming chromosomes from Roman numerals to Arabic
raw_oriC['Organism'] = raw_oriC['Organism'].str.replace('chromosome IV', 'chromosome 4', regex = True)
raw_oriC['Organism'] = raw_oriC['Organism'].str.replace('chromosome III', 'chromosome 3', regex = True)
raw_oriC['Organism'] = raw_oriC['Organism'].str.replace('chromosome II', 'chromosome 2', regex = True)
raw_oriC['Organism'] = raw_oriC['Organism'].str.replace('chromosome I', 'chromosome 1', regex = True)

#####################################################################################################################################
### Katelyn found that organisms with multiple oriC usually flank the DnaA region and wondered whether they were
### really two separate oriC or just one that was broken up by the gene. I'll make two dataframes. One that keeps
### it the way it is now, with each multi-oriC organism getting its own entry, and one that merges two (or more)
### oriC if they flank the same DnaA region.
#####################################################################################################################################

# First dataset: everything is kept like it was
oriC_sep = raw_oriC.copy() # sep = separate
oriC_sep.to_csv('DoriC_oriC_split.csv', index=False)

# Second dataset: merge DnaA gene if flanked by two oriC
multi_oriC = raw_oriC[raw_oriC['RefSeq'].duplicated(keep=False) == True].copy()
del raw_oriC

dups_refseq = multi_oriC.pivot_table(columns=['RefSeq'], aggfunc='size')
more_than_two_oriC = dups_refseq[dups_refseq > 2]
# print(f'There are {len(more_than_two_oriC.values.tolist())} organsims for which DoriC predicted more than two oriC')
# print(f'I will delete these, because there are only 5 and the merging with flanking DnaA makes my brain hurt.')

two_oriC = multi_oriC.copy()
del multi_oriC
for i in more_than_two_oriC.keys().to_list():
    two_oriC = two_oriC[two_oriC['RefSeq'] != i]

# Changing oriC_f to DnaA_f if oriC_f flows into DnaA_i 
two_oriC.reset_index(inplace=True)
filtered_two_oriC = two_oriC.copy()
filtered_two_oriC
for sample in two_oriC.iterrows():
    i = sample[0]
    df = sample[1] # sample is a tuple with (index from older version of , values)
    if i > 0 and two_oriC.iloc[i-1]['RefSeq'] == df['RefSeq']:
        # Option 1: oriC_1 -> DnaA -> oriC_2
        if df['DnaA_0_i'] - 1 == df['oriC_f'] and two_oriC.iloc[i-1]['oriC_i'] == two_oriC.iloc[i-1]['DnaA_0_f'] + 1:
            filtered_two_oriC.loc[filtered_two_oriC['RefSeq'] == df['RefSeq'], 'oriC_i'] = df['oriC_i']
            filtered_two_oriC.loc[filtered_two_oriC['RefSeq'] == df['RefSeq'], 'oriC_f'] = two_oriC.iloc[i-1]['oriC_f']
        # Option 2: oriC_2 -> DnaA -> oriC_1
        elif df['DnaA_0_f'] + 1 == df['oriC_i'] and two_oriC.iloc[i-1]['oriC_f'] == two_oriC.iloc[i-1]['DnaA_0_i'] - 1:
            filtered_two_oriC.loc[filtered_two_oriC['RefSeq'] == df['RefSeq'], 'oriC_i'] = two_oriC.iloc[i-1]['oriC_i']
            filtered_two_oriC.loc[filtered_two_oriC['RefSeq'] == df['RefSeq'], 'oriC_f'] = df['oriC_f']

# This change does not merge two oriCs if the DnaA starts at position 1 or ends at the last position of the sequence
# Because the total length of the sequence is not known, it is not possible to know if an oriC ends at the last position in the sequence.
# For an example check: NZ_CP024969
filtered_two_oriC.drop('OriC sequence', axis=1, inplace=True) # Sequence is not right anymore.
filtered_two_oriC = filtered_two_oriC.drop_duplicates(subset=['RefSeq'])

# Can continue work on this once on the cluster.