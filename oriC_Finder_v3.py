# Libraries
import scipy.signal as sp
import numpy as np
import pandas as pd
import os
from itertools import combinations
from itertools import product

# Self-made module
import plotting_functions as pf

# Set cwd to location of this script
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )


def read_FASTA(filename):
    """Read in a file in FASTA format and returns the name and sequence"""
    seq = ''
    name = None
    with open(filename, 'r') as fh:
        for line in fh:
            line = line.rstrip()
            if len(line) == 0 and name is not None:
                break # In case of empty line at end of file
            elif line[0] == '>' and name is None:
                name = line[1:]
            elif line[0] == '>' and name is not None:
                break # In case the FASTA contains multiple sequences, only read the first one.
            else:
                seq = seq + line
    return name, seq.upper()


def calc_everything(seq):
    """
    Z-curve, GC-skew, and N-count calculation. In one function so only one iteration of the sequence is necessary.
    Input:
        seq     : string DNA sequence
    Return:
        x, y, z : 1D-np.arrays of the three Z-curve components
        gc      : 1D-np.array of the GC-skew
        n       : N-count of the sequence. N can be any nucleotide.
    """
    x, y, z = [], [], []
    gc_skew = []
    a, c, t, g, n = 0, 0, 0, 0, 0

    for base in seq:
        if base == "A":
            a +=1                  
        elif base == "C":
            c +=1
        elif base == "G":
            g +=1
        elif base == "T":
            t +=1 
        elif base == "N":
            n += 1

        gc_skew.append(g - c)
        x.append( (a + g)-(c + t) ) # Purine vs Pyrimidine
        y.append( (a + c)-(g + t) ) # Amino vs Keto
        z.append( (a + t)-(c + g) ) # Weak vs Strong Hydrogen Bonds
    return np.asarray(x), np.asarray(y), np.asarray(z), np.asarray(gc_skew), n


def detect_peaks(curve):
    """Calculates peaks of 1D-np.array and returns its indeces."""
    maxima, _ = sp.find_peaks( curve, distance=len(curve)//12)
    maxima    = np.append(maxima, curve.argmax())
    minima, _ = sp.find_peaks( np.negative(curve), distance=len(curve)//12)
    minima    = np.append(minima, curve.argmin())
    return np.unique(np.concatenate( (maxima, minima), axis=0))


def get_peak_windows(curve_size, peaks, window_size=500):
    """Calculates the indeces of the windows around each peak of a curve with given size"""
    windows = []
    frame   = window_size//2

    for i in range(len(peaks)):
        down = peaks[i] + frame
        up   = peaks[i] - frame
        # This is for cases where the extremes are at the beginning and end of the domain
        if up < 0:
            a_up = [j for j in range(0, peaks[i])]
            b_up = [j for j in range(curve_size-1 + up, curve_size)]
        else:
            a_up = [j for j in range(up, peaks[i]+1)]
            b_up = None
        
        if down > curve_size-1:
            a_down = [j for j in range(peaks[i], curve_size)]
            b_down = [j for j in range(0, down - curve_size)]
        else:
            a_down = [j for j in range(peaks[i], down+1)]
            b_down = None

        up_window = a_up + b_up if b_up is not None else a_up
        down_window = a_down + b_down if b_down is not None else a_down
        up_window.extend(down_window)
        window = sorted(up_window)
        windows.append( window )
    return windows


def split_window(window):
    '''Splits a peak_window in two based on whether the indeces in the window are consecutive'''
    in_a = True
    a, b = [], []
    for i, val in enumerate(window):
        if i-1 > 0 and val != window[i-1] + 1:
            in_a = False
        a.append(val) if in_a else b.append(val)
    return a, b


def filter_peaks(curve, peaks, peak_windows, mode='max'):
    """
    Filters the given peaks based on the type of extreme it should be and the area around the peak in two steps.
        Filter 1: Check if any windows intersect the window of another peak (includes circularity of DNA)
        Filter 2: Check if peaks are actually the extreme in their windows
    Input:
        curve          : 1D-np.array
        peaks          : peaks of curve
        peak_windows   : peak_windows of peaks
        mode           : 'max'|'min'. Which type of extreme do you want to find?
    Return:
        accepted_peaks : peaks that passed both filters
    """
    rejected_peaks = []
    for(i, win_i), (j, win_j) in combinations(enumerate(peak_windows), 2):
        # Filter 1: Check if any other windows intersecting the window of peak i
        if win_i != win_j and curve[peaks[i]] != curve[peaks[j]] and len(set(win_i).intersection(win_j)) > 0:
            # Either peaks at the beginning and end of the DNA sequence or a peak found by sp.find_peaks() that is very close to the global min/max 
            if mode == 'max':
                rejected_peaks.extend( np.where(curve == min(curve[peaks[i]], curve[peaks[j]]) )[0].tolist() )
            elif mode == 'min':
                rejected_peaks.extend( np.where(curve == max(curve[peaks[i]], curve[peaks[j]]) )[0].tolist() )

    for i, win_i in enumerate(peak_windows):
        # Filter 2: Check if peaks are actually the extreme in their windows
        if len(curve)-1 in win_i and 0 in win_i:
            a, b = split_window(win_i)
            comparator_win = a if peaks[i] in a else b
        else:
            comparator_win = win_i
    
        if mode == 'max' and np.max( curve[comparator_win] ) > curve[peaks[i]]:
                rejected_peaks.append(peaks[i])
        elif mode == 'min' and np.min( curve[comparator_win] ) < curve[peaks[i]]:
                rejected_peaks.append(peaks[i])

    # Create list of peaks that passed both filters
    rejected_peaks = list(set(rejected_peaks))
    accepted_peaks = [x for x in peaks if x not in rejected_peaks]

    # Fail-safe in case no peaks were accepted, I don't expect this to ever be needed, since global extremes never get rejected
    if len(accepted_peaks) == 0:
        print('it was needed... that\'s bad...')
        if mode == 'max':
            accepted_peaks.append(np.argmax(curve))
        elif mode == 'min':
            accepted_peaks.append(np.argmin(curve))
    return accepted_peaks


def get_peaks_to_merge(peaks, peak_windows):
    """
    Get the indeces that give the same value in the curve and are in eachothers window.
    Input:
        peaks          : peaks of a curve
        peak_windows   : peak_windows of peaks
    Return:
        peaks_to_merge : list of tuples. Each tuple contains the two indeces that have to be merged
    """
    peaks_to_merge = []
    for i, win_i in enumerate(peak_windows):
        for j in range(len(peaks)):
            if peaks[j] in win_i and i != j:
                # No need to check if they give the same value on the curve, because if they didn't, one would have been rejected by filter_peaks() already
                peaks_to_merge.append( tuple( sorted([peaks[i], peaks[j]]) ) )
    return list(set(peaks_to_merge))


def merge_peaks(curve_size, a, b):
    """
    Calculate the distance between index a and b on circular DNA
    Returns the middle point between the shortest distance
    """
    dist_1 = max(a, b) - min(a, b)
    dist_2 = min(a, b) + curve_size-1 - max(a, b)
    if dist_2 < dist_1:
        merged_idx = min(a, b) - dist_2//2
        if merged_idx < 0:
            merged_idx = curve_size-1 + merged_idx
    else:
        merged_idx = min(a, b) + dist_1//2
    return merged_idx


def match_peaks(peaks_x, peaks_y, peak_windows_x, peak_windows_y):
    """
    Checks if the peaks from x line up with the peaks from y.
    If they don't, then consult the gc-skew.
    """
    matched_peaks = []
    for (i_x, win_x), (j_y, win_y) in product(enumerate(peak_windows_x), enumerate(peak_windows_y)):
        if len( set(win_x).intersection(win_y) ) > 0:
            matched_peaks.append( (peaks_x[i_x], peaks_y[j_y]) )
    return list(set(matched_peaks))


def process_array(curve, mode='max', window_size=500):
    """
    Runs the given 1D-array (curve) through all processing functions for oriC identification.
    Returns its peaks and windows.
    """
    init_peaks = detect_peaks(curve)
    init_peak_windows = get_peak_windows(curve.size, init_peaks, window_size=window_size)
    accept_peaks = filter_peaks(curve, init_peaks, init_peak_windows, mode=mode)
    accept_windows = get_peak_windows(curve.size, accept_peaks, window_size=window_size)
    peaks_to_merge = get_peaks_to_merge(accept_peaks, accept_windows)

    # Necessary, because global extreme can be a plateau and there might be two peaks only 5 bp apart
    single_peaks = [x for x in accept_peaks if not any(x in y for y in peaks_to_merge)]
    merged_peaks = [merge_peaks(curve.size, to_merge[0], to_merge[1]) for to_merge in peaks_to_merge]
    peaks = single_peaks + merged_peaks
    peak_windows = get_peak_windows(curve.size, peaks, window_size=window_size)
    return peaks, peak_windows


def get_false_order(seq, curve, oriC_locations, mode='max', window_size=500):
    '''
    Get n_penalty score for how an oriC was predicted. The higher the penalty, the lower the prediction's reliability.
        - N-count in oriC           : if high and multiple oriC, could be a false extreme 
        - N-count in whole sequence : if high, could mean that the oriC has not yet been sequenced
    Return:
        false_order : T|F, whether the order of the oriCs could be wrong due to 'N' bases in area around oriC
    '''
    false_orders = []
    if len(oriC_locations) > 1:
        windows = get_peak_windows(len(seq), oriC_locations, window_size=window_size)
        n_per_oriC = []
        for window in windows:
            n = 0
            for i in window:
                n += 1 if seq[i] == 'N' else 0
            n_per_oriC.append(n)

        # NOTE: only checks the oriC in top position against the other oriCs, not every combination
        # This looks at cases where there are two peaks of similar heights/depths that are quite far apart.
        # e.g. if the x-curve was W-shaped, making it so you have two very similar minima. (or y-curve like M)
        for i in range(1, len(n_per_oriC)):
            if mode == 'max':
                false_orders.append(curve[oriC_locations[0]] - n_per_oriC[0] <= curve[oriC_locations[i]] + n_per_oriC[i])
            elif mode == 'min':
                false_orders.append(curve[oriC_locations[0]] + n_per_oriC[0] >= curve[oriC_locations[i]] - n_per_oriC[i])
    return any(false_orders)


def calc_dist(curve_size, a, b):
    # same distance calculation as in merge_peaks()
    dist_1 = max(a, b) - min(a, b)
    dist_2 = min(a, b) + curve_size-1 - max(a, b)
    return min(dist_1, dist_2)


def curve_combinations(curves_list, peaks_list, windows_list):
    '''Get every matched_peaks combination for x, y, and gc.'''
    oriC_locations_list = []
    for (i, peaks_i), (j, peaks_j) in combinations(enumerate(peaks_list), 2):
        matched_peaks  = match_peaks(peaks_i, peaks_j, windows_list[i], windows_list[j])
        oriC_locations_list.append( [merge_peaks(len(curves_list[0]), matches[0], matches[1]) for matches in matched_peaks] )
    return oriC_locations_list


def get_adj_mat(curve_size, peaks):
    '''Gets adjacency matrix for given peaks'''
    adj_mat = np.zeros((len(peaks), len(peaks)))
    for (i_a, a), (i_b, b) in combinations(enumerate(peaks), r=2):
        dist = calc_dist(curve_size, a, b)
        adj_mat[i_a, i_b] = dist
        adj_mat[i_b, i_a] = dist
    return adj_mat


def get_connected_groups(peaks, adj_mat, threshold):
    visited = [False] * len(peaks)
    connected_groups_idx = []
    for i in range(len(peaks)):
        if not visited[i]:
            group = []
            _, _, visited, group, _ = _DFS_recurse(i, adj_mat, visited, group, threshold=threshold)
            connected_groups_idx.append(group)
    connected_groups_vals = [ [peaks[i] for i in idx_group] for idx_group in connected_groups_idx ]
    return connected_groups_vals


def _DFS_recurse(idx, adj_mat, visited, connected_list, threshold):
    visited[idx] = True
    connected_list.append(idx)
    for i in range(len(visited)):
        if i == idx:
            continue
        elif adj_mat[i][idx] <= threshold and not visited[i]:
            _, _, visited, connected_list, _ = _DFS_recurse(i,adj_mat,visited, connected_list, threshold)
    return idx, adj_mat, visited, connected_list, threshold


def merge_oriCs(curve_size, groups):
    '''Finds the average index of a group and returns those values. groups is a nested-list'''
    mutable = sorted( groups, key=lambda x:len(x), reverse=True )
    oriCs, occurances = [], []
    for group in mutable:
        group.sort()
        for i in range(len(group)):
            if group[-1] - group[i] >= (curve_size-1)/2:
                group[i] += curve_size-1

        avg_val = sum(group)//len(group)
        if avg_val > curve_size-1:
            avg_val -= curve_size-1
        oriCs.append(avg_val)
        occurances.append(len(group))

    total_pot_oriCs = len( [y for x in mutable for y in x] )
    occurances = [x/total_pot_oriCs for x in occurances]
    return oriCs, occurances


def find_oriCs(filename):
    '''
    VERSION 3
    Locates potential oriC based on Z-curve and GC-skew analysis.
    Three window_sizes are used: 1, 3 and 5 % of the total genome length. The oriCs that were found by most combinations, get returned.
    Input:
        filename    : FASTA-file with circular bacterial genome or chromosome
        oriC_size   : The size of the predicted oriCs
        window_size : Affects the minimal distance between potential oriC
    Return:
        properties  : Dictionary with oriC properties
    '''
    # NOTE: Potential further improvements:
    #   - AT% measure of oriC: should be characteristically higher than the rest of the genome.
    #   - Essential gene proximity: genes like DnaA are very likely to be close to or on the oriC
    name, sequence = read_FASTA(filename)
    x, y, z, gc, n = calc_everything(sequence)

    windows = [0.01, 0.03, 0.05]
    peaks   = []

    for fraction in windows:
        window_size = int(len(sequence) * fraction)
        peaks_x,  peak_windows_x  = process_array(x , mode='min', window_size=window_size)
        peaks_y,  peak_windows_y  = process_array(y , mode='max', window_size=window_size)
        peaks_gc, peak_windows_gc = process_array(gc, mode='min', window_size=window_size)
        peaks.extend( [y for x in curve_combinations( (x, y, gc), (peaks_x, peaks_y, peaks_gc), (peak_windows_x, peak_windows_y, peak_windows_gc) ) for y in x] )

    # Connected components in undirected graph problem
    matrix = get_adj_mat(len(sequence), peaks)

    # Depth-First Search
    connected_groups = get_connected_groups(peaks, matrix, int(len(sequence)*windows[-1]))

    # Remove potential oriC if it was not matched to any other.
    connected_groups = [x for x in connected_groups if len(x) > 1]

    # Get oriCs
    oriCs, occurances = merge_oriCs(len(sequence), connected_groups)

    # Check false order
    false_order = any( get_false_order( sequence, i[0], oriCs, mode=i[1], window_size=int(len(sequence)*windows[-1]) ) for i in ((x, 'min'), (y, 'max'), (gc, 'min')) )        

    # Final dictionary
    preferred_oriC_properties = {
        'name'         : name,
        'oriC_middles' : oriCs,
        'occurances'   : occurances,
        'z_curve'      : (x, y, z),
        'gc_skew'      : gc,
        'false_order'  : false_order,
        'seq_size'     : len(sequence),
        'gc_conc'      : ( sequence.count('G') + sequence.count('C') ) / len(sequence)
    }

    return preferred_oriC_properties


if __name__ == '__main__':
    # For testing a small folder
    # for fasta in os.listdir('./test_fastas'):
    #     file = os.path.join('test_fastas', fasta)
    #     preferred_properties = find_oriCs(file)

    #     name    = preferred_properties['name']
    #     Z_curve = preferred_properties['z_curve']
    #     GC_skew = preferred_properties['gc_skew']

    #     print(name)
    #     print('QoP  :', preferred_properties['occurances'])
    #     print('oriCs:', preferred_properties['oriC_middles'])

    #     pf.plot_Z_curve_2D(list(Z_curve[:2]) + [GC_skew], [preferred_properties['oriC_middles']]*3, name)
    #     # pf.plot_skew(GC_skew, [preferred_properties['oriC_middles']], name)
    #     # pf.plot_Z_curve_3D(Z_curve, name)

    # For Testing single files
    properties = find_oriCs('./worst_fastas/v3_min_0.6/NZ_CP018845.fasta')
    name    = properties['name']
    Z_curve = properties['z_curve']
    GC_skew = properties['gc_skew']

    print(name)
    print('QoP  :', properties['occurances'], properties['false_order'])
    print('oriCs:', properties['oriC_middles'])


    pf.plot_Z_curve_2D(list(Z_curve[:2]) + [GC_skew], [properties['oriC_middles']]*3, name)
    # pf.plot_skew(GC_skew, [properties['oriC_middles']], name)
    # pf.plot_Z_curve_3D(Z_curve, name)
