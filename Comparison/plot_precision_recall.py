import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PRECISION   = 0
RECALL      = 0
DISTANCE    = 0
PREC_VS_REC = 1

df = pd.read_csv('Comparison\All_occurances_precision_recall.csv')

x = df['min_confidence']
y_names = ['Precision (%)', 'Recall (%)', 'Distance off (%)']

y_precision = df['precision_A'], df['precision_Z'], df['precision_G'], df['precision_D']
y_recall    = df['recall_A'], df['recall_Z'], df['recall_G'], df['recall_D']
y_distance  = df['distance_pc_A'], df['distance_pc_Z'], df['distance_pc_G'], df['distance_pc_D']

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.88)

if PRECISION:
    p1_1, = ax.plot(x, y_precision[0], label='Avg')
    p1_2, = ax.plot(x, y_precision[1], label='Z')
    p1_3, = ax.plot(x, y_precision[2], label='G')
    p1_4, = ax.plot(x, y_precision[3], label='D')
    ax.set_ylabel(y_names[0])
    ax.set_title('Precision graph of ZoriC on DoriC dataset.\nTP when accurate within 2.5% of genome length')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 105, 5)])
    ax.set_ylim(0, 100)


if RECALL:
    p1_1, = ax.plot(x, y_recall[0], label='Avg')
    p1_2, = ax.plot(x, y_recall[1], label='Z')
    p1_3, = ax.plot(x, y_recall[2], label='G')
    p1_4, = ax.plot(x, y_recall[3], label='D')
    ax.set_ylabel(y_names[1])
    ax.set_title('Recall graph of ZoriC on DoriC dataset.\nTP when accurate within 2.5% of genome length')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 105, 5)])
    ax.set_ylim(0, 100)

if DISTANCE:
    p1_1, = ax.plot(x, y_distance[0], label='Avg')
    p1_2, = ax.plot(x, y_distance[1], label='Z')
    p1_3, = ax.plot(x, y_distance[2], label='G')
    p1_4, = ax.plot(x, y_distance[3], label='D')
    ax.set_ylabel(y_names[2])
    ax.set_title('Average distance off graph of ZoriC on DoriC dataset.')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 21)])
    ax.set_ylim(0, 20)

if PREC_VS_REC:
    p1_1 = ax.scatter(y_recall[0], y_precision[0], label='Avg')
    p1_2 = ax.scatter([86.52], [83.95], c='r', marker='x')
    p1_3 = ax.plot([0, 45], [83.95, 83.95], 'r--')
    p1_3 = ax.plot([45, 45], [0, 83.95], 'r--')
    p1_3 = ax.plot([0, 86.52], [44, 44], 'r-.')
    p1_3 = ax.plot([86.52, 86.52], [0, 44], 'r-.')
    ax.set_yticks([x for x in range(0, 110, 10)])
    ax.set_ylim(0, 100)
    ax.set_xlabel("Recall (%)")
    ax.set_ylabel("Precision (%)")

ax.set_xlim(0, 100)
tkw = dict(size=4, width=1.5, gridOn=True)
ax.tick_params(axis='y', **tkw)
ax.tick_params(axis='x', **tkw)
if not PREC_VS_REC:
    ax.yaxis.label.set_color(p1_1.get_color())
    ax.set_xlabel("Minimum confidence (%)")
    ax.legend(handles=[p1_1, p1_2, p1_3, p1_4])

plt.show()