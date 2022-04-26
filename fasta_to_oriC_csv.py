import os
# Set cwd to location of this script
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

from oriC_Finder import find_oriCs

import pandas as pd
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 30)

path = os.path.join( os.getcwd(), 'refseq_15', 'bacteria')
sample_list = os.listdir( path )

df_dict = {
    'RefSeq'  : [],
    'Organism': [],
    'oriC_i'  : [],
    'oriC_f'  : [],
    'QoP'     : [] # Quality of Prediction
}

for fasta in sample_list:
    raw_name, oriC_ranges, _, _, QoP = find_oriCs( os.path.join(path, fasta) )
    
    # Name processing
    name_list = raw_name.split(' ')
    RefSeq = name_list[0]
    Organism = ' '.join(name_list[1:])

    # OriC processing
    for i, oriC in enumerate(oriC_ranges):

    # Appending to dataframe
        df_dict['RefSeq'].append(RefSeq)
        df_dict['Organism'].append(Organism)
        df_dict['oriC_i'].append(oriC[0])
        df_dict['oriC_f'].append(oriC[1])
        df_dict['QoP'].append(QoP)

df = pd.DataFrame.from_dict(df_dict)

# Removing version numbers for now
df['RefSeq'] = df['RefSeq'].str.extract(r'([^.]*)')

print(df.head(10))

df.to_csv(f'NCBI_oriC_{df.shape[0]}_improved.csv', index=False)