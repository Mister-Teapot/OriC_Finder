import pandas as pd
import matplotlib.pyplot as plt

# Pandas printing options
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 20)

CSV_1 = 'v1'
CSV_2 = 'v3'

distances_v1 = pd.read_csv('Comparison/'+CSV_1+'/hist_all_df.csv').sort_values(by=['RefSeq'], ignore_index=True)
distances_load = pd.read_csv('Comparison/'+CSV_2+'/hist_all_df.csv').sort_values(by=['RefSeq'], ignore_index=True)

if distances_load.shape[0] > distances_v1.shape[0]:
    distances_v2 = distances_load[distances_load['RefSeq'].isin(distances_v1['RefSeq'])].copy()
    distances_v2.reset_index(drop=True, inplace=True)
    distances_v2.sort_values(by=['RefSeq'], ignore_index=True, inplace=True)
else:
    distances_v2 = distances_load.copy()
    distances_v1 = distances_v1[distances_v1['RefSeq'].isin(distances_v2['RefSeq'])].copy()
    distances_v1.reset_index(drop=True, inplace=True)
    distances_v1.sort_values(by=['RefSeq'], ignore_index=True, inplace=True)
del distances_load

df = pd.read_csv('Comparison/v2/in_both_sets_all.csv').sort_values(by=['RefSeq'], ignore_index=True)

# best i could come up with
for i, sample in distances_v2.iterrows():
    if i == 0:
        df_w_dups = df[df['RefSeq'] == sample['RefSeq']]
    else:
        df_w_dups = pd.concat([df_w_dups, df[df['RefSeq'] == sample['RefSeq']]], axis=0, ignore_index=True)

merged = pd.concat([df_w_dups, distances_v2], axis=1)
merged = merged.T.drop_duplicates().T.copy()

better_distance = distances_v2[distances_v2['Distance'] < distances_v1['Distance'].values]
same_distance = distances_v2[distances_v2['Distance'] == distances_v1['Distance'].values]
worse_distance = distances_v2[distances_v2['Distance'] > distances_v1['Distance'].values]

rel_better = distances_v1[distances_v1.index.isin(better_distance.index)]['Distance'] - better_distance['Distance']
rel_worse = worse_distance['Distance'] - distances_v1[distances_v1.index.isin(worse_distance.index)]['Distance']

rel_better.rename('Difference', inplace=True)
rel_worse.rename('Difference', inplace=True)

sorted_better = pd.concat([better_distance, rel_better], axis=1)
sorted_worse = pd.concat([worse_distance, rel_worse], axis=1).sort_values(by=['Difference'])
# print(sorted_worse.tail())

print('All organisms                  :', distances_v2.shape[0])
print('Organisms with better distance :', better_distance.shape[0])
print('Organisms with same distance   :', same_distance.shape[0])
print('Organisms with worse distance  :', worse_distance.shape[0])
print()
print('Distances improved on average :', int(rel_better.mean()), 'bp')
print('Distances worsened on average :', int(rel_worse.mean()), 'bp')
print()
print('Average seq_len of genomes with improved predictions :', int(merged[merged['RefSeq'].isin(better_distance['RefSeq'])]['Sequence_length'].mean()), 'bp')
print('Average seq_len of genomes with worsened predictions :', int(merged[merged['RefSeq'].isin(worse_distance['RefSeq'])]['Sequence_length'].mean()), 'bp')
# print()
# print(f'Average option of genomes with improved predictions : {merged[merged["RefSeq"].isin(better_distance["RefSeq"])]["O_penalty"].mean():.3f}')
# print(f'Average option of genomes with worsened predictions : {merged[merged["RefSeq"].isin(worse_distance["RefSeq"])]["O_penalty"].mean():.3f}')

plt.grid(axis='both', which='both')
plt.scatter(merged['Distance'], merged['Sequence_length'])
plt.xlabel('Distance from oriC (bp)')
plt.ylabel('Sequence length (bp)')
plt.show()

plt.grid(axis='both', which='both')
plt.scatter( merged['Distance']/merged['Sequence_length']*100, merged['Sequence_length'])
plt.xlabel('Distance from oriC (% of total genome)')
plt.ylabel('Sequence length (bp)')
plt.show()