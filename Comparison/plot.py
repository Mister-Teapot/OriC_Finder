import pandas as pd
import matplotlib.pyplot as plt

summary_both = pd.read_csv('Comparison/v3/Summary.csv')
# summary_only = pd.read_csv('Comparison/v3/Summary_only_false_orders.csv')
# summary_none = pd.read_csv('Comparison/v3/Summary_no_false_order.csv')

names = ['DoriCs found (%)', 'Correct ZoriCs (%)', 'Distance off (%)', 'Distance off (bp)', '# Organisms for which at least one ZoriC was found']
x = summary_both['Confidence']

y_1 = summary_both['DoriC oriCs found'], summary_both['Z_oriCs correct'], summary_both['Distance_pc'], summary_both['Distance_bp'], summary_both['accession_with_a_ZoriC']
# y_2 = summary_none['DoriC oriCs found'], summary_none['Z_oriCs that found a DoriC oriC'], summary_none['Distance_percentage'], summary_none['Distance_bp'], summary_none['Organisms with a ZoriC']
# y_3 = summary_only['DoriC oriCs found'], summary_only['Z_oriCs that found a DoriC oriC'], summary_only['Distance_percentage'], summary_only['Distance_bp'], summary_only['Organisms with a ZoriC']

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.75)

p1_1, = ax.plot(x, y_1[0], "b-", label='All Samples')
# p1_2, = ax.plot(x, y_2[2], "g--", label='No False order')
# p1_3, = ax.plot(x, y_3[2], "g-.", label='Only False order')

twin1 = ax.twinx()
p2_1, = twin1.plot(x, y_1[1], "r-", label='All Samples')
# p2_2, = twin1.plot(x, y_2[1], "r--", label='No False order')
# p2_3, = twin1.plot(x, y_3[1], "r-.", label='Only False order')

twin2 = ax.twinx()
twin2.spines.right.set_position(("axes", 1.2))
p3_1, = twin2.plot(x, y_1[2], "g-", label='All Samples')
# p3_2, = twin2.plot(x, y_2[2], "g--", label='No False order')
# p3_3, = twin2.plot(x, y_3[2], "g-.", label='Only False order')

ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
twin1.set_ylim(0, 100)
twin2.set_ylim(0, 10)

ax.set_xlabel("Minimum confidence (%)")
ax.set_ylabel(names[0])
twin1.set_ylabel(names[1])
twin2.set_ylabel(names[2])

ax.yaxis.label.set_color(p1_1.get_color())
twin1.yaxis.label.set_color(p2_1.get_color())
twin2.yaxis.label.set_color(p3_1.get_color())

tkw = dict(size=4, width=1.5)
ax.tick_params(axis='y', colors=p1_1.get_color(), **tkw)
twin1.tick_params(axis='y', colors=p2_1.get_color(), **tkw)
twin2.tick_params(axis='y', colors=p3_1.get_color(), **tkw)
ax.tick_params(axis='x', **tkw)

# ax.legend(handles=[p1_1, p1_2, p1_3])

plt.show()



'''
v1 = pd.read_csv(r'C:\0. School\Bachelor Thesis\Zoya_Code+Data\GitHub_Repositories\OriC_Finder\Comparison\v1\hist_all_df.csv').sort_values(by=['Distance'], ignore_index=True)
v2 = pd.read_csv(r'C:\0. School\Bachelor Thesis\Zoya_Code+Data\GitHub_Repositories\OriC_Finder\Comparison\v2\hist_all_df.csv').sort_values(by=['Distance'], ignore_index=True)
v3 = pd.read_csv(r'C:\0. School\Bachelor Thesis\Zoya_Code+Data\GitHub_Repositories\OriC_Finder\Comparison\v3\hist_all_df.csv').sort_values(by=['Distance'], ignore_index=True)

best_1  = v1.head(20)['RefSeq']
worst_1 = v1.tail(20)['RefSeq']

best_2  = v2.head(20)['RefSeq']
worst_2 = v2.tail(20)['RefSeq']

best_3  = v3.head(20)['RefSeq']
worst_3 = v3.tail(20)['RefSeq']

print('Overlap best/worst: v1 -> v2:')
print(f'{len(set(best_1).intersection(best_2)) + len(set(worst_1).intersection(worst_2))} out of {len(worst_2)*2} are the same.')
print('Overlap best/worst: v1 -> v3:')
print(f'{len(set(best_1).intersection(best_3)) + len(set(worst_1).intersection(worst_3))} out of {len(worst_3)*2} are the same.')
print('Overlap best/worst: v2 -> v3:')
print(f'{len(set(best_2).intersection(best_3)) + len(set(worst_3).intersection(worst_3))} out of {len(worst_3)*2} are the same.')
'''