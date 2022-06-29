import pandas as pd
import matplotlib.pyplot as plt

PRECISION = 0
RECALL    = 1
DISTANCE  = 0

df = pd.read_csv('Comparison/version_3_4_5_precision_recall.csv')

x = df['min_confidence']
y_names = ['Precision (%)', 'Recall (%)', 'Distance off (%)']

y_precision = df['precision_v3'], df['precision_v4'], df['precision_v5']
y_recall    = df['recall_v3'], df['recall_v4'], df['recall_v5']
y_distance  = df['distance_pc_v3'], df['distance_pc_v4'], df['distance_pc_v5']

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.88)

if PRECISION:
    p1_1, = ax.plot(x, y_precision[0], "b--", label='V3')
    p1_2, = ax.plot(x, y_precision[1], "b-", label='V4')
    p1_3, = ax.plot(x, y_precision[2], "b-.", label='V5')
    ax.set_ylabel(y_names[0])
    ax.set_title('Precision graph of ZoriC on DoriC dataset.\nTP when accurate within 2.5% of genome length')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(60, 105, 5)])
    ax.set_ylim(60, 100)


if RECALL:
    p1_1, = ax.plot(x, y_recall[0], "r--", label='V3')
    p1_2, = ax.plot(x, y_recall[1], "r-", label='V4')
    p1_3, = ax.plot(x, y_recall[2], "r-.", label='V5')
    ax.set_ylabel(y_names[1])
    ax.set_title('Recall graph of ZoriC on DoriC dataset.\nTP when accurate within 2.5% of genome length')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 105, 5)])
    ax.set_ylim(0, 100)

if DISTANCE:
    p1_1, = ax.plot(x, y_distance[0], "g--", label='V3')
    p1_2, = ax.plot(x, y_distance[1], "g-", label='V4')
    p1_3, = ax.plot(x, y_distance[2], "g-.", label='V5')
    ax.set_ylabel(y_names[2])
    ax.set_title('Average distance off graph of ZoriC on DoriC dataset.')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 11)])
    ax.set_ylim(0, 10)


ax.set_xlim(0, 100)
ax.set_xlabel("Minimum confidence (%)")
ax.yaxis.label.set_color(p1_1.get_color())
tkw = dict(size=4, width=1.5, gridOn=True)
ax.tick_params(axis='y', colors=p1_1.get_color(), **tkw)
ax.tick_params(axis='x', **tkw)
ax.legend(handles=[p1_1, p1_2, p1_3])

plt.show()