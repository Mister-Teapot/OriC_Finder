import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PRECISION   = 0
RECALL      = 0
DISTANCE    = 0
PREC_VS_REC = 0

df_no = pd.read_csv('Comparison\All_occurances_precision_recall_no_G.csv')
df_ml = pd.read_csv('Comparison\precision_recall_standard_no_G_sep_G.csv')

x_no = df_no['min_confidence']
x_ml = df_ml['min_confidence']
y_names = ['Precision (%)', 'Recall (%)', 'Distance off (%)']

y_precision_no = df_no['precision_A'], df_no['precision_Z'], df_no['precision_G'], df_no['precision_D']
y_recall_no    = df_no['recall_A'], df_no['recall_Z'], df_no['recall_G'], df_no['recall_D']
y_distance_no  = df_no['distance_pc_A'], df_no['distance_pc_Z'], df_no['distance_pc_G'], df_no['distance_pc_D']

y_precision_ml = df_ml['precision_standard'], df_ml['precision_no_G'], df_ml['precision_sep_G']
y_recall_ml    = df_ml['recall_standard'], df_ml['recall_no_G'], df_ml['recall_sep_G']
y_distance_ml  = df_ml['distance_pc_standard'], df_ml['distance_pc_no_G'], df_ml['distance_pc_sep_G']


fig, ax = plt.subplots()
tw = ax.twiny()
fig.subplots_adjust(right=0.88, top=0.8)

if PRECISION:
    p1_1, = ax.plot(x_no, y_precision_no[3], 'r', label='Mean')
    p1_2, = tw.plot(x_ml, y_precision_ml[2], 'b', label='SVC')
    ax.set_ylabel(y_names[0])
    ax.set_title('Precision graph of ORCA on DoriC dataset.\nTP when accurate within 2.5% of genome length')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(30, 105, 5)])
    ax.set_ylim(30, 100)


if RECALL:
    p1_1, = ax.plot(x_no, y_recall_no[3], 'r', label='Mean')
    p1_2, = tw.plot(x_ml, y_recall_ml[2], 'b', label='SVC')
    ax.set_ylabel(y_names[1])
    ax.set_title('Recall graph of ORCA on DoriC dataset.\nTP when accurate within 2.5% of genome length')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 120, 10)])
    ax.set_ylim(0, 110)

if DISTANCE:
    p1_1, = ax.plot(x_no, y_distance_no[3], 'r', label='Mean')
    p1_2, = tw.plot(x_ml, y_distance_ml[2], 'b', label='SVC')
    ax.set_ylabel(y_names[2])
    ax.set_title('Average distance off graph of ORCA on DoriC dataset.')
    ax.set_xticks([x for x in range(0, 110, 10)])
    ax.set_yticks([x for x in range(0, 18, 2)])
    ax.set_ylim(0, 16)

if PREC_VS_REC:
    p1_1 = ax.scatter(y_recall_no[3], y_precision_no[3], c='r', label='Mean')
    p1_2 = ax.scatter(y_recall_ml[2], y_precision_ml[2], c='b', label='SVC')
    ax.set_yticks([x for x in range(0, 110, 10)])
    ax.set_ylim(0, 110)
    ax.set_xlabel("Recall (%)", fontsize=15)
    ax.set_ylabel("Precision (%)", fontsize=15)

ax.set_xlim(0, 100)
tw.set_xlim(-2, 3)
tw.set_xticks([x/10 for x in range(-20, 35, 5)])
tkw = dict(size=4, width=1.5, gridOn=True)
ax.tick_params(axis='y', **tkw)
ax.tick_params(axis='x', **tkw)
if not PREC_VS_REC:
    ax.set_xlabel("Minimum confidence (%)")
    tw.set_xlabel("Minimum confidence")
    ax.xaxis.label.set_color(p1_1.get_color())
    tw.xaxis.label.set_color(p1_2.get_color())
    ax.legend(handles=[p1_1, p1_2])

plt.show()