import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


palette = sns.color_palette(["r", "b"], 2)
palette2 = sns.color_palette(["maroon", "darkblue"], 2)
column = 'Prediction'

# DATASET   = [# Accessions that have been experimentally verified.
#     'NC_000964', 'NC_002947', 'NC_003272', 'NC_003869', 'NC_005090', 'NC_006461', 'NC_007633', 'NC_000913', 'NC_003047',
#     'NC_007604', 'NC_000962', 'NC_002696', 'NC_002971', 'NC_005363', 'NC_008255', 'NC_009850', 'NC_010546', 'NC_011916'
# ]

# df = pd.read_csv('Hyperparameter tuning/tuning_hist_plot_SVC.csv')
df = pd.read_csv('Hyperparameter tuning/tuning_hist_plot_exp.csv')
accs = []
for i, acc in df['RefSeq_oriC'].iteritems():
    accs.append(acc[:-2])
print(len(set(accs)))
counts = 0
for i in set(accs):
    counts += accs.count(i)
print(counts/len(set(accs)))

# df[column] = (df[column] - df[column].min()) / (df[column].max() - df[column].min())   
# df = pd.concat([df, new_versions], axis=0).reset_index()
ax = sns.histplot(x=df[column], hue=df['Correct'], element='step', palette=palette, fill=True, bins=25, **{'linewidth': 2})
# ax2 = sns.histplot(x=new_versions[column], hue=new_versions['Correct'], element='step', palette=palette2, fill=True, bins=25, **{'linewidth': 2})
# ax = sns.histplot(x=df['Avg_occurance'], hue=df['Correct'], element='step', multiple='layer', palette=palette, fill=True, bins=100, **{'linewidth': 2})

TP = df.loc[(df['Correct'] == True) & (df[column] >= 0)].shape[0]
TN = df.loc[(df['Correct'] == False) & (df[column] < 0)].shape[0]
FP = df.loc[(df['Correct'] == False) & (df[column] >= 0)].shape[0]
FN = df.loc[(df['Correct'] == True) & (df[column] < 0)].shape[0]

print(f'precision : {TP / (TP+FP)*100:.2f}')
print(f'recall    : {TP / (TP+FN)*100:.2f}')
print('TP:', TP)
print('FN:', FN)

legend = ax.get_legend()
# legend2 = ax2.get_legend()
handles = legend.legendHandles
# handles2 = legend2.legendHandles
legend.remove()
# legend2.remove()
leg = ax.legend(handles, ['False', 'True'], loc='upper left')
# leg2 = ax2.legend(handles2, ['F', 'T'], loc='upper right')
ax.add_artist(leg)
# ax2.add_artist(leg2)

plt.ylabel('True Positive/Negative Count', fontsize=15)
plt.xlabel('Confidence', fontsize=15)
plt.yticks([x for x in range(0, 30, 2)])
# plt.grid(True)
plt.ylim(0, 28)
plt.xlim(-3, 3)
plt.show()