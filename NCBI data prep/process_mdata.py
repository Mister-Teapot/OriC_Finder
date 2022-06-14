import pandas as pd
pd.set_option('display.width', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)

table = pd.read_table('NCBI data prep/mdata.tab')

for i, sample in table.iterrows():
    if isinstance(sample['Accession'], float):
        acc = ''
    else:
        acc = sample['Accession'].split(' ')[0]
    table.at[i, 'Accession'] = acc
table['Accession'] = table['Accession'].str.extract(r'([^.]*)')
table.to_csv('NCBI data prep/mdata_clean.csv', index=False)

refseqs = table[table['RefSeq'] == True]['Accession']
for i in refseqs:
    print(i)