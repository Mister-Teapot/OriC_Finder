import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

shift = 4000000
df = pd.read_csv('shift.csv').to_numpy().T[1:,:]

fig = plt.figure(figsize=(10,4))
fig.subplots_adjust(bottom=0.29, top=0.82, left=0.08, right=0.97)
ax = plt.axes()
tw = ax.twiny()

ax.set_xlabel('Sequence index (bp)')
ax.grid(True, which='major', zorder=1)

ax.plot([len(df[0]), len(df[0])], [df[0][-1], df[0][0]], c='r', linestyle='--', zorder=2)
tw.plot([0, 0], [df[1][-1], df[1][0]], c='b', linestyle='--', zorder=2)

ax.plot(range(len(df[0])), df[0], c='r', label='$y_n$ (0 bp shift)', zorder=3)
ax.plot(range(len(df[0]), len(df[0])+shift), df[0][:shift], c='r', zorder=3)

tw.plot(range(len(df[1])), df[1], c='b', label='$y_n$ (4 Mbp shift)', zorder=3)
tw.plot(range(-shift, 0), df[1][-shift:], c='b', zorder=3)


ax.scatter(
    [np.argmax(df[0]), len(df[0])],
    [np.max(df[0]), np.max(df[0])],
    c='darkred',marker='$\mathrm{O}$', s=60, label='Maximum (0 bp shift)', zorder=4
)
tw.scatter(
    [len(df[1])-shift, -shift],
    [df[1][len(df[1])-shift], df[1][-shift]],
    c='darkred',marker='$\mathrm{O}$', s=60, zorder=4
)
tw.scatter(
    [np.argmax(df[1]), np.argmax(df[1])],
    [np.max(df[1]), df[0][np.argmax(df[1])+shift]],
    c='darkblue', marker='$\mathrm{O}$', s=60, label='Maximum (4 Mbp shift)', zorder=4
)

ax.set_xlim(0, shift+len(df[1]))
tw.set_xlim(-shift, len(df[1]))

ax.ticklabel_format(axis='x', style='sci', scilimits=(3,3), useMathText=True)
tw.ticklabel_format(axis='x', style='sci', scilimits=(3,3), useMathText=True)

ax.tick_params(axis ='x', colors='r')
tw.tick_params(axis ='x', colors='b')

h1, l1 = ax.get_legend_handles_labels()
h2, l2 = tw.get_legend_handles_labels()

plt.legend(
    handles=h1+h2,
    labels=l1+l2,
    bbox_to_anchor=(0.2, -0.45, 0.5, .102),
    loc='center',
    ncol=2,
    mode="expand",
    borderaxespad=0.
)
plt.show()