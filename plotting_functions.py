import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_Z_curve_3D(Z_curve, name):
    """3D-plot function with name as tilte"""
    fig = plt.figure(figsize=(6,6))
    ax = plt.axes(projection='3d')
    ax.set_xlabel('Purine vs. Pyrimidine', fontsize=10)
    ax.set_ylabel('Amino vs. Keto', fontsize=10)
    ax.set_zlabel('Weak vs. Strong H-bond', fontsize=10)

    x, y, z = Z_curve
    ax.plot3D(x, y, z, linewidth=0.8)
    ax.set_title(f'Z-Curve: {name}', fontsize=10,loc='center', pad=20)
    plt.show()


def plot_Z_curve_2D(y_val_list, peaks, name):
    """
    Plots 2D Z-curve. Can display up to 4 y-axes in a single figure.
        y_val_list : list of lists with y-axis values
        peaks      : list of lists with indeces of peaks for arrays in y_val_list
        name       : used in plot title
    """
    color_list = ['r', 'b', 'g', 'k']
    fig = plt.figure(figsize=(7,5), tight_layout=True)
    base_ax = plt.axes()
    ax_list = [base_ax] + [base_ax.twinx() for i in range(len(y_val_list) - 1)]

    offset = 1
    for axis in ax_list[1:]:
        axis.spines.right.set_position(("axes", offset))
        offset += 0.2

    for i, ax in enumerate(y_val_list):
        peaks_y = np.asarray([ax[j] for j in peaks[i]]) # y refers to the y-axis coordinates, not the y-curve
        ax_list[i].plot(range(len(ax)), ax, color_list[i], zorder=1)
        ax_list[i].scatter(peaks[i], peaks_y, marker='X', c='k', zorder=2)
        ax_list[i].tick_params(axis ='y', colors=color_list[i])

    base_ax.set_title(f'2D Z-Curve: {name}', fontsize=10,loc='center', pad=20)
    plt.show()


def plot_skew(skewArray, peaks, name):
    """Plots single skew diagram and its peaks"""

    fig = plt.figure()
    ax1 = plt.axes()

    peaks_y = np.asarray([skewArray[i] for i in peaks])

    ax1.set_title(f'GC-skew: {name}', fontsize=10,loc='center', pad=20)
    ax1.plot(range(len(skewArray)), skewArray, 'r', zorder=1)
    ax1.scatter(peaks, peaks_y, marker='X', c='k', zorder=2)
    plt.show()


def distance_histogram(db, log=False):
    plt.hist(db['Distance_bp'], bins=[x for x in range(0, 3600000//2, 1000)], log=log) # 3.6e6 = avg. len of bacterial chromosome
    plt.show()
    plt.hist(db['Distance_pc'], bins=[x for x in range(0, 50, 1)], log=log)
    plt.show()