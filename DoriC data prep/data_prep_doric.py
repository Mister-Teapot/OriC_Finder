import re
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

Roughly speaking, we're comparing the Z-curve method against the GC-skew method
(DoriC does not purely use the Z-curve but also combines it with GC-skew)
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

# Remove RefSeq versions in accension number (might keep them later, since DoriC is not up to date with NCBI in 2022)
raw_oriC['RefSeq'] = raw_oriC['RefSeq'].str.extract(r'([^.]*)')
raw_oriC['Lineage'] = raw_oriC['Lineage'].str.rstrip('.')

# NOT removing multichromosomal organisms -> Organisms with multiple chromosomes already have separate entries

# Renaming chromosomes from Roman numerals to Arabic
raw_oriC['Organism'] = raw_oriC['Organism'].str.replace('chromosome IV', 'chromosome 4', regex = True)
raw_oriC['Organism'] = raw_oriC['Organism'].str.replace('chromosome III', 'chromosome 3', regex = True)
raw_oriC['Organism'] = raw_oriC['Organism'].str.replace('chromosome II', 'chromosome 2', regex = True)
raw_oriC['Organism'] = raw_oriC['Organism'].str.replace('chromosome I', 'chromosome 1', regex = True)


#####################################################################################################################################

### The GC-skew method will only predict 1 oriC. It will be interesting to see if the oriC it predicts aligns
### with ANY oriC found by the DoriC team. We want to try to incorporate multi-oriC organisms anyway, so this
### is a start ¯\_(ツ)_/¯
### Katelyn found that organisms with multiple oriC usually flank the DnaA region and wondered whether they were
### really two separate oriC or just one that was broken up by the gene. I'll make two dataframes. One that keeps
### it the way it is now, with each multi-oriC organism getting its own entry, and one that merges two (or more)
### oriC if they flank the same DnaA region.

oriC_sep = raw_oriC.copy() # sep = separate
oriC_sep.to_csv('DoriC_oriC_split.csv', index=False)

###################################################################################################################################
oriC_lengths_DoriC = [ i[1]['oriC_f'] - i[1]['oriC_i'] for i in oriC_sep.iterrows() if i[1]['oriC_f'] - i[1]['oriC_i'] > 0 ]

oriC_lengths_eperiments = []
with open(r'C:\0. School\Bachelor Thesis\ref_lengths.txt', 'r') as fh:
    for line in fh:
        line = line.strip()
        i, f = line.split('..')
        oriC_lengths_eperiments.append(float(f) - float(i))

print(f'Average length of oriC in DoriC                 : {sum(oriC_lengths_DoriC)/len(oriC_lengths_DoriC):.2f}')
print(f'Average length of oriC experimentally           : {sum(oriC_lengths_eperiments)/len(oriC_lengths_eperiments):.2f}')
print(f'Not taking into account oriC\'s that run 3\' -> 5\': {(1 - len(oriC_lengths_DoriC)/oriC_sep.shape[0]) * 100:.2f} %')
###################################################################################################################################

# duplicates = raw_oriC[raw_oriC['RefSeq'].duplicated(keep=False) == True]
