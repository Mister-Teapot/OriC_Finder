import pandas as pd
from Bio import Entrez
import multiprocessing as mp
pd.set_option('display.width', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)

table = pd.read_table('NCBI data prep/mdata.tab')
refseqs = []
for i, sample in table.iterrows():
    if isinstance(sample['ChromosomeAccession'], str):
        accs = [ acc.split('.')[0] for acc in sample['ChromosomeAccession'].split(' ') if '_' in acc]
        refseqs.extend(accs)

Entrez.email = 'zoyavanmeel@gmail.com'
Entrez.api_key = '795d705fb638507c9b2295c89cc64ee88108'

test = pd.read_table('NCBI data prep/test.tab')
a = set(table['Assembly'].to_list()).intersection(test)
print(len(a), a)

# def do(i):
#     print(f'Processing:\t{refseqs[i]}\t{i+1}/{len(refseqs)}')
#     handle = Entrez.efetch(db='nuccore', id=refseqs[i], rettype='gb_with_parts', retmode='text')
#     for j in handle:
#         if 'rep_' in j or 'REP_' in j:
#             print('Found rep_origin in:', refseqs[i])

# if __name__ == '__main__':
#     with mp.Pool(mp.cpu_count() - 1) as pool:
#         pool.map(do, [x for x in range(len(refseqs))])
# # print(len(refseqs))