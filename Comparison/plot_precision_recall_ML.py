import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PRECISION   = 1
RECALL      = 0
DISTANCE    = 0
PREC_VS_REC = 0

df = pd.read_csv('Comparison\precision_recall_standard_no_G_sep_G.csv')

x = df['min_confidence']
y_names = ['Precision (%)', 'Recall (%)', 'Distance off (%)']

y_precision = df['precision_standard'], df['precision_no_G'], df['precision_sep_G']
y_recall    = df['recall_standard'], df['recall_no_G'], df['recall_sep_G']
y_distance  = df['distance_pc_standard'], df['distance_pc_no_G'], df['distance_pc_sep_G']

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.88)

if PRECISION:
    p1_1, = ax.plot(x, y_precision[0], c='r', linestyle='-', linewidth=2, label='Standard')
    p1_2, = ax.plot(x, y_precision[1], c='b', linestyle='--', linewidth=2, label='No-G')
    p1_3, = ax.plot(x, y_precision[2], c='g', linestyle='-.', linewidth=2, label='Separate-G')
    # p1_4, = ax.plot([0, 0], [0, 100], 'r', label='decision boundary')
    ax.set_ylabel(y_names[0], fontsize=15)
    # ax.set_title('Precision graph of ZoriC on DoriC dataset.\nTP when accurate within 2.5% of geNo-e length')
    low, high, step = 0, 110, 10
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(low, high+step, step)])
    ax.set_ylim(low, high)


if RECALL:
    p1_1, = ax.plot(x, y_recall[0], c='r', linestyle='-', linewidth=2, label='Standard')
    p1_2, = ax.plot(x, y_recall[1], c='b', linestyle='--', linewidth=2, label='No-G')
    p1_3, = ax.plot(x, y_recall[2], c='g', linestyle='-.', linewidth=2, label='Separate-G')
    # p1_4, = ax.plot([0, 0], [0, 100], 'r', label='decision boundary')
    ax.set_ylabel(y_names[1], fontsize=15)
    # ax.set_title('Recall graph of ZoriC on DoriC dataset.\nTP when accurate within 2.5% of geNo-e length')
    low, high, step = 0, 110, 10
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(low, high+step, step)])
    ax.set_ylim(low, high)

if DISTANCE:
    p1_1, = ax.plot(x, y_distance[0], c='r', linestyle='-', linewidth=2, label='Standard')
    p1_2, = ax.plot(x, y_distance[1], c='b', linestyle='--', linewidth=2, label='No-G')
    p1_3, = ax.plot(x, y_distance[2], c='g', linestyle='-.', linewidth=2, label='Separate-G')
    # p1_4, = ax.plot([0, 0], [0, 100], 'r', label='decision boundary')
    ax.set_ylabel(y_names[2], fontsize=15)
    # ax.set_title('Average distance off graph of ZoriC on DoriC dataset.')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 21)])
    ax.set_ylim(0, 5)

if PREC_VS_REC:
    p1_1 = ax.scatter(y_recall[0], y_precision[0], label='Standard')
    p1_1 = ax.scatter(y_recall[1], y_precision[1], label='No-G')
    p1_1 = ax.scatter(y_recall[2], y_precision[2], label='Separate-G')
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
    ax.set_xlim(-2, 3)
    ax.set_xticks([x/10 for x in range(-20, 35, 5)])
    ax.set_xlabel("Minimum confidence", fontsize=15)
    ax.legend(handles=[p1_1, p1_2, p1_3])

plt.show()