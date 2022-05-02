import os, sys
# Set cwd to location of this script
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )
sys.path.append('../OriC-Finder')

from OriC-Finder import oriC_Finder 
import pandas as pd
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 30)

path = os.path.join( os.getcwd(), 'refseq_15', 'bacteria')
sample_list = os.listdir( path )

df_dict = {
    'RefSeq'         : [],
    'Organism'       : [],
    'Sequence_length': [],
    'Penalties'      : [],
    'False_order'    : []
}

for fasta in sample_list:
    properties = oriC_Finder.find_oriCs( os.path.join(path, fasta) )
    
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
    for i in range(len(properties['oriC_middles'])):
        df_dict.update( {f'oriC_edges_{i}'  : properties['oriC_edges'][i]  } )
        df_dict.update( {f'oriC_middles_{i}': properties['oriC_middles'][i]} )

df = pd.DataFrame.from_dict(df_dict)

# Removing version numbers for now
df['RefSeq'] = df['RefSeq'].str.extract(r'([^.]*)')

print(df.head(10))

# df.to_csv(f'NCBI_oriC_{df.shape[0]}_improved.csv', index=False)