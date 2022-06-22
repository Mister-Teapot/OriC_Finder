import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('Comparison/version_3_and_4_precision_recall.csv')

x = df['min_confidence']
y_names = ['Precision (%)', 'Recall (%)']

y_precision = df['precision_v3'], df['precision_v4']
y_recall    = df['recall_v3'], df['recall_v4']

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.88)

# p1_1, = ax.plot(x, y_precision[0], "b--", label='V3')
# p1_2, = ax.plot(x, y_precision[1], "b-", label='V4')

p1_1, = ax.plot(x, y_recall[0], "r--", label='V3')
p1_2, = ax.plot(x, y_recall[1], "r-", label='V4')

ax.set_xlim(0, 100)
ax.set_xticks([x for x in range(0, 110, 10)])
ax.set_ylim(40, 100)
ax.set_yticks([x for x in range(40, 110, 5)])

ax.set_xlabel("Minimum confidence (%)")
# ax.set_ylabel(y_names[0])
ax.set_ylabel(y_names[1])

ax.yaxis.label.set_color(p1_1.get_color())

tkw = dict(size=4, width=1.5, gridOn=True)
ax.tick_params(axis='y', colors=p1_1.get_color(), **tkw)
ax.tick_params(axis='x', **tkw)
ax.set_title('Recall graph of ZoriC on DoriC dataset.\nTP when accurate within 2.5% of genome length')
ax.legend(handles=[p1_1, p1_2])

plt.show()