import pandas as pd

df_3k = pd.read_csv(r'C:\0. School\Bachelor Thesis\Zoya_Code+Data\OriFinder\OriC_Finder\NCBI data prep\refseq_8\NCBI_oriC_7.csv')
df_rest = pd.read_csv(r'C:\0. School\Bachelor Thesis\Zoya_Code+Data\OriFinder\OriC_Finder\NCBI data prep\3k+299+15_merged.csv')
df_3k.dropna(how='all', axis=1, inplace=True)
df_rest.dropna(how='all', axis=1, inplace=True)
df_3k['RefSeq'] = df_3k['RefSeq'].str.extract(r'([^.]*)')

df_merged = pd.concat((df_3k, df_rest), axis=0, ignore_index=True)

df_merged.to_csv(r'C:\0. School\Bachelor Thesis\Zoya_Code+Data\OriFinder\OriC_Finder\NCBI data prep\3k+299+15+7_merged.csv', index=False)