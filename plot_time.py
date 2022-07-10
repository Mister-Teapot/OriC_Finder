import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import seaborn as sns

sys.path.append('../OriC_Finder/')

df = pd.read_csv('Time_performance.csv').sort_values(by='seq_len', ascending=True)
df['Index'] = df['Index'].sort_values()
df['Index'] = df['Index'].astype(str)

sns.set_theme(style="whitegrid")
f, ax = plt.subplots(figsize=(12, 4))

# ax.bar(df["seq_len"], df["total_time"], color='#4455ff', width=80000)
# ax.bar(df["seq_len"], df['calc_disp_time']+df['read_genes_time'], color='#7582ff', width=80000)
# ax.bar(df["seq_len"], df["calc_disp_time"], color='#d6d9ff', width=80000)

# ax.step(df["seq_len"], df["total_time"])

# ax.set_xlim(0, 14000000)
sns.set_color_codes("bright")
sns.barplot(y="total_time", x="seq_len", data=df,
    label="Score calculation", color="b")

sns.set_color_codes("pastel")
sns.barplot(y=df['calc_disp_time']+df['read_genes_time'], x=df['seq_len'],
    label="Gene file parsing", color="b")

sns.set_color_codes("muted")
sns.barplot(y="calc_disp_time", x="seq_len", data=df,
    label="Disparity calculation", color="b")


mean = df['total_time'].mean()
print(df['total_time'].mean())
print(df['calc_disp_time'].mean())
print(df['read_genes_time'].mean())
print(df['seq_len'].mean())
print(df['num_of_genes'].mean())
ax.plot([0, 102], [mean, mean], c='b')
ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5,1.15,0,0), frameon=True)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize=6)
ax.set(xlim=(-.5, 101.5),
    xlabel="",
    ylabel="Prediction time (s)",
    yticks=[x for x in range(0, 110, 10)],
    ylim=(0,100)
)
sns.despine(left=True, bottom=True)
# ax = sns.displot(data=df, x='seq_len', palette='b', bins=25)
# plt.ticklabel_format(axis='x', style='sci', useMathText=True, scilimits=(6,6))
plt.show()

# fig, ax = plt.subplots()
# fig.subplots_adjust(right=0.825)

# ax.set_xlabel('Sequence length (bp)', fontsize=15)
# ax.set_ylabel('Genes', fontsize=15)

# ax.ticklabel_format(axis='y', style='sci', useMathText=True)

# ax.set_xlim(0, 14000000)
# ax.set_ylim(0, 10000)
# ax.set_xticks([x for x in range(0, 15400000, 1400000)])
# ax.set_yticks([x for x in range(0, 11000, 1000)])

# x = df['seq_len'].to_numpy()
# y = df['num_of_genes'].to_numpy()

# slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
# a = ax.scatter(df['seq_len'], df['num_of_genes'], c='r', marker='x', label='Original data')
# print('a:', slope*10000, 'b:', intercept)
# b = ax.plot(x, slope*x + intercept, 'b', label='$y = 8.32 \cdot 10^{-4} \cdot x + 311.48$')
# print('R**2', r_value**2)
# ax.grid(True)

# ax.legend()
# plt.show()
