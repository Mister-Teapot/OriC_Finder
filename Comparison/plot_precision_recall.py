import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PRECISION   = 0
RECALL      = 0
DISTANCE    = 0
PREC_VS_REC = 1

df = pd.read_csv('Comparison\All_occurances_precision_recall.csv')
df_2 = pd.read_csv('Comparison\All_occurances_precision_recall_no_G.csv')

x = df['min_confidence']
y_names = ['Precision (%)', 'Recall (%)', 'Distance off (%)']

y_precision = df['precision_A'], df['precision_Z'], df['precision_G'], df['precision_D'], df_2['precision_A']
y_recall    = df['recall_A'], df['recall_Z'], df['recall_G'], df['recall_D'], df_2['recall_A']
y_distance  = df['distance_pc_A'], df['distance_pc_Z'], df['distance_pc_G'], df['distance_pc_D'], df_2['distance_pc_A']

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.88)

if PRECISION:
    p1_1, = ax.plot(x, y_precision[0], c='r', linestyle='-', linewidth=2, label='Standard')
    p1_2, = ax.plot(x, y_precision[4], c='b', linestyle='--', linewidth=2, label='No-G')
    # p1_2, = ax.plot(x, y_precision[1], c='b', linestyle='--', linewidth=2, label='Z-score')
    # p1_3, = ax.plot(x, y_precision[2], c='g', linestyle='-.', linewidth=2, label='G-score')
    # p1_4, = ax.plot(x, y_precision[3], c='k', linestyle=':', linewidth=2, label='D-score')
    ax.set_ylabel(y_names[0], fontsize=15)
    # ax.set_title('Precision graph of ZoriC on DoriC dataset.\nTP when accurate within 2.5% of genome length')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 110, 10)])
    ax.set_ylim(0, 100)


if RECALL:
    p1_1, = ax.plot(x, y_recall[0], c='r', linestyle='-', linewidth=2, label='Standard')
    p1_2, = ax.plot(x, y_recall[4], c='b', linestyle='--', linewidth=2, label='No-G')
    # p1_2, = ax.plot(x, y_recall[1], c='b', linestyle='--', linewidth=2, label='Z-score')
    # p1_3, = ax.plot(x, y_recall[2], c='g', linestyle='-.', linewidth=2, label='G-score')
    # p1_4, = ax.plot(x, y_recall[3], c='k', linestyle=':', linewidth=2, label='D-score')
    ax.set_ylabel(y_names[1], fontsize=15)
    # ax.set_title('Recall graph of ZoriC on DoriC dataset.\nTP when accurate within 2.5% of genome length')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 110, 10)])
    ax.set_ylim(0, 100)

if DISTANCE:
    p1_1, = ax.plot(x, y_distance[0], c='r', label='Mean')
    p1_2, = ax.plot(x, y_distance[1], c='b', label='Z-score')
    p1_3, = ax.plot(x, y_distance[2], c='g', label='G-score')
    p1_4, = ax.plot(x, y_distance[3], c='k', label='D-score')
    ax.set_ylabel(y_names[2], fontsize=15)
    # ax.set_title('Average distance off graph of ZoriC on DoriC dataset.')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 21)])
    ax.set_ylim(0, 20)

if PREC_VS_REC:
    p1_1 = ax.scatter(y_recall[4], y_precision[4], c='b', marker='^', label='Mean-ORCA')
    p1_2 = ax.scatter([90.05], [93.64], c='r', marker='x', label='SVC-ORCA')
    p1_3 = ax.plot([90.05, 90.05], [0, 93.64], 'r--')
    p1_3 = ax.plot([0, 90.05], [93.64, 93.64], 'r--')
    ax.set_yticks([x for x in range(0, 110, 10)])
    ax.set_ylim(30, 100)
    ax.set_xlabel("Recall (%)", fontsize=15)
    ax.set_ylabel("Precision (%)", fontsize=15)
    ax.legend(handles=[p1_1, p1_2]) #, p1_3, p1_4


ax.set_xlim(30, 100)
tkw = dict(size=4, width=1.5, gridOn=True)
ax.tick_params(axis='y', **tkw)
ax.tick_params(axis='x', **tkw)
if not PREC_VS_REC:
    ax.set_xlabel("Minimum confidence (%)", fontsize=15)
    ax.legend(handles=[p1_1, p1_2]) #, p1_3, p1_4

plt.show()