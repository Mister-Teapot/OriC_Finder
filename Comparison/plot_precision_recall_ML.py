import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PRECISION   = 1
RECALL      = 0
DISTANCE    = 0
PREC_VS_REC = 0

df = pd.read_csv('Comparison\precision_recall_standard_no_G.csv')

x = df['min_confidence']
y_names = ['Precision (%)', 'Recall (%)', 'Distance off (%)']

y_precision = df['precision_standard'], df['precision_no_G']
y_recall    = df['recall_standard'], df['recall_no_G']
y_distance  = df['distance_pc_standard'], df['distance_pc_no_G']

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.88)

if PRECISION:
    p1_1, = ax.plot(x, y_precision[0], label='standard')
    p1_2, = ax.plot(x, y_precision[1], label='no G')
    ax.set_ylabel(y_names[0])
    ax.set_title('Precision graph of ZoriC on DoriC dataset.\nTP when accurate within 2.5% of genome length')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 105, 5)])
    ax.set_ylim(0, 100)


if RECALL:
    p1_1, = ax.plot(x, y_recall[0], label='standard')
    p1_2, = ax.plot(x, y_recall[1], label='no G')
    ax.set_ylabel(y_names[1])
    ax.set_title('Recall graph of ZoriC on DoriC dataset.\nTP when accurate within 2.5% of genome length')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 105, 5)])
    ax.set_ylim(0, 100)

if DISTANCE:
    p1_1, = ax.plot(x, y_distance[0], label='standard')
    p1_2, = ax.plot(x, y_distance[1], label='no G')
    ax.set_ylabel(y_names[2])
    ax.set_title('Average distance off graph of ZoriC on DoriC dataset.')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 21)])
    ax.set_ylim(0, 5)

if PREC_VS_REC:
    p1_1 = ax.scatter(y_recall[0], y_precision[0], label='standard')
    p1_1 = ax.scatter(y_recall[1], y_precision[1], label='no G')
    ax.set_yticks([x for x in range(0, 110, 10)])
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_ylim(0, 100)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Recall (%)")
    ax.set_ylabel("Precision (%)")

tkw = dict(size=4, width=1.5, gridOn=True)
ax.tick_params(axis='y', **tkw)
ax.tick_params(axis='x', **tkw)
if not PREC_VS_REC:
    ax.set_xlim(-1.5, 1)
    ax.set_xticks([x/10 for x in range(-15, 15, 5)])
    ax.set_xlabel("Minimum confidence (%)")
    ax.legend(handles=[p1_1, p1_2])

plt.show()