# Libraries
import scipy.signal as sp
import numpy as np
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


def detect_peaks(skew_array):
    """Calculates peaks of 1D-np.array and returns its indeces."""
    maxima, _ = sp.find_peaks( skew_array, distance=len(skew_array)//12)
    maxima    = np.append(maxima, skew_array.argmax())
    minima, _ = sp.find_peaks( np.negative(skew_array), distance=len(skew_array)//12)
    minima    = np.append(minima, skew_array.argmin())
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
    If that still does not give anything: consult Jasmijn
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


def process_matches(curves_list, peaks_list, windows_list):

    # Option 0: Base oriC position on x and y components of the Z-curve
    matched_peaks = match_peaks(peaks_list[0], peaks_list[1], windows_list[0], windows_list[1])
    option = 0

    if len(matched_peaks) == 0:
        # Option 1: Base oriC position solely on global extremes with a larger window_size
        matched_peaks = get_last_resort_positions(curves_list[0], curves_list[1], window_size=len(windows_list[0][0])*2)
        option += 1

    if len(matched_peaks) == 0:
        # Option 2: Base oriC position of peaks of x and gc-skew
        matched_peaks = match_peaks(peaks_list[0], peaks_list[2], windows_list[0], windows_list[2])
        option += 1

    if len(matched_peaks) == 0:
        # Option 3: Base oriC position of peaks of y and gc-skew
        matched_peaks = match_peaks(peaks_list[1], peaks_list[2], windows_list[1], windows_list[2])
        option += 1

    if len(matched_peaks) == 0:
        # Option 4: Base oriC position of global minimum of x
        matched_peaks = [(curves_list[0].argmin(), curves_list[0].argmin())]
        option += 1

    return matched_peaks, option


def get_last_resort_positions(x, y, window_size=500):
    """Get min of x, max of y and their respective peak_windows"""
    x_peak, y_peak = [x.argmin()], [y.argmax()]
    x_window = get_peak_windows(len(x), x_peak, window_size=window_size)
    y_window = get_peak_windows(len(y), y_peak, window_size=window_size)
    return match_peaks(x_peak, y_peak, x_window, y_window)


def sort_oriCs(curve_a, curve_b, oriC_locations, mode_a='max', mode_b='max', window_size=500):
    """
    Sorts the found oriC locations based on prominence in the Z-curve and GC-skew.
    Best oriC is (close to) global min/max of given curves
    The rest is sorted by the size of the codomain if the curves with the peak_window of the location as its domain.
    Essentially, if peak around oriC is bigger, it is better.
    Does not sort in-place
    """
    extreme_a = curve_a.argmax() if mode_a == 'max' else curve_a.argmin()
    extreme_b = curve_b.argmax() if mode_b == 'max' else curve_b.argmin()
    prominences = []
    for i, oriC_pos in enumerate(oriC_locations):
        oriC_window = get_peak_windows(curve_a.size, [oriC_pos], window_size=window_size)[0]
        if extreme_a in oriC_window or extreme_b in oriC_window:
            prominences.append( (i, float('inf')) )
        else:
            # Calculate average peak prominence in both curves
            a = np.max(curve_a[oriC_window]) - np.min(curve_a[oriC_window])
            b = np.max(curve_b[oriC_window]) - np.min(curve_b[oriC_window])
            avg = (a + b)/2
            prominences.append( (i, avg) )
    prominences.sort(reverse=True, key=lambda x: x[1])
    return [oriC_locations[x[0]] for x in prominences]


def get_oriC_ranges(seq_len, oriC_locations, range_size=500):
    """Get the start and end index of the given locations"""
    # window = list of all indeces in a range
    # range  = tuple of two edges extracted from a window
    full_oriC_windows = get_peak_windows(seq_len, oriC_locations, window_size=range_size)
    oriC_ranges = []
    for oriC_window in full_oriC_windows:
        if seq_len-1 in oriC_window and 0 in oriC_window:
            # oriC on domain borders -> split into a and b to work with.
            a, b = split_window(oriC_window)
            # oriC_ranges.append( (min(a), max(a)) )
            # oriC_ranges.append( (min(b), max(b)) )
            oriC_ranges.append( (min(b), max(a)) )
        else:
            oriC_ranges.append( (oriC_window[0], oriC_window[-1]) )
    return oriC_ranges


def get_false_order(seq, curve_a, curve_b, oriC_locations, mode_a='max', mode_b='max', window_size=500):
    '''
    Get n_penalty score for how an oriC was predicted. The higher the penalty, the lower the prediction's reliability.
        - N-count in oriC           : if high and multiple oriC, could be a false extreme 
        - N-count in whole sequence : if high, could mean that the oriC has not yet been sequenced
    Return:
        n_penalties : n_penalty score
        false_order : T|F, whether the order of the oriCs could be wrong due to 'N' bases in area around oriC
    '''
    false_order = False
    if len(oriC_locations) > 1:
        windows = get_peak_windows(len(seq), oriC_locations, window_size=window_size)
        n_per_oriC = []
        for window in windows:
            n = 0
            for i in window:
                n += 1 if seq[i] == 'N' else 0
            n_per_oriC.append(n)

        # NOTE: only checks the oriC in global extreme against the other oriCs, not every combination
        # This looks at cases where there are two peaks of similar heights/depths that are quite far apart.
        # e.g. if the x-curve was W-shaped, making it so you have two very similar minima. (or y-curve like M)
        false_orders = [0]
        for i in range(1, len(n_per_oriC)):
            pot_false_order = 0 # potential false ordering: not reliable
            if mode_a == 'max' and curve_a[oriC_locations[0]] - n_per_oriC[0] <= curve_a[oriC_locations[i]] + n_per_oriC[i]:
                pot_false_order += 1
            elif mode_a == 'min' and curve_a[oriC_locations[0]] + n_per_oriC[0] >= curve_a[oriC_locations[i]] - n_per_oriC[i]:
                pot_false_order += 1
            if mode_b == 'max' and curve_b[oriC_locations[0]] - n_per_oriC[0] <= curve_b[oriC_locations[i]] + n_per_oriC[i]:
                pot_false_order += 1
            elif mode_b == 'min' and curve_b[oriC_locations[0]] + n_per_oriC[0] >= curve_b[oriC_locations[i]] - n_per_oriC[i]:
                pot_false_order += 1
            false_orders.append(pot_false_order)
        # Filter out non-false orders
        false_orders = [x for x in false_orders if x != 0]
        if len(false_orders) > 0: # This is bad: unreliable result
            false_order = True
    return false_order


def get_dist_penalty(curve_size, matched_peaks, dist_penalty_factor=10):
    # same distance calculation as in merge_peaks()
    penalties = []
    for match in matched_peaks:
        dist_1 = max(match[0], match[1]) - min(match[0], match[1])
        dist_2 = min(match[0], match[1]) + curve_size-1 - max(match[0], match[1])
        penalties.append(min(dist_1, dist_2))
    return [i * dist_penalty_factor for i in penalties]


def curve_combinations(curves_list, peaks_list, windows_list):
    '''Get every matched_peaks combination for x, y, and gc.'''
    oriC_locations_list = []
    for (i, peaks_i), (j, peaks_j) in combinations(enumerate(peaks_list), 2):
        mode_a = 'max' if i == 1 else 'min' # min for x and gc, max for y
        mode_b = 'max' if j == 1 else 'min'

        matched_peaks  = match_peaks(peaks_i, peaks_j, windows_list[i], windows_list[j])
        oriC_locations = [merge_peaks(len(curves_list[0]), matches[0], matches[1]) for matches in matched_peaks]
        oriC_locations_list.append( sort_oriCs(curves_list[i], curves_list[j], oriC_locations, mode_a=mode_a, mode_b=mode_b, window_size=len(windows_list[0][0])) )
    return oriC_locations_list

# On second thought: Make this a processing and oriC finding function for two given curves. Then you can use this function in find_oriCs() and import this function elsewhere too.
# TODO: Make this function and get all options for all data to determine which options are best-suited when.
# TODO: Also check the distance calculation (NC_016946 inspection)
# TODO: Also check when data is too noisy to be properly interpreted.


def find_oriCs(filename, oriC_size=500):#, window_size=150000):
    '''
    Locates potential oriC based on Z-curve and GC-skew analysis.
    Input:
        filename    : FASTA-file with circular bacterial genome or chromosome
        oriC_size   : The size of the predicted oriCs
        window_size : Affects the minimal distance between potential oriC
    Return:
        name        : Name of the organism. Read from FASTA-file.
        oriC_ranges : List of tuples with the indeces of the potential oriC in the given FASTA-file. 
                      tuple = (first_idx, last_idx). The oriC are ordered on importance.
                      Importance is based on the severity of the polarity change around the oriC peak.
        Z_Curve     : Tuple of the x, y, z components of the Z-curve analysis
        GC_skew     : Values of the GC-skew analysis
        option      : The option that was used to determine the oriCs. 0 = best; 4 = worst.
    '''
    name, sequence = read_FASTA(filename)
    x, y, z, gc, n = calc_everything(sequence)
    false_order    = False
    window_size    = int(len(sequence)*0.03)

    # z gives no information for oriC location prediction
    peaks_x,  peak_windows_x  = process_array(x , mode='min', window_size=window_size)
    peaks_y,  peak_windows_y  = process_array(y , mode='max', window_size=window_size)
    peaks_gc, peak_windows_gc = process_array(gc, mode='min', window_size=window_size)
    pf.plot_Z_curve_2D([x, y], [peaks_x, peaks_y], name)

    # Getting prefered option
    matched_peaks, option = process_matches( (x, y, gc), (peaks_x, peaks_y, peaks_gc), (peak_windows_x, peak_windows_y, peak_windows_gc) )
    oriC_locations = [merge_peaks(len(sequence), matches[0], matches[1]) for matches in matched_peaks]

    # Sort oriC locations based on importance and calculate penalties
    if option == 0:
        oriC_locations = sort_oriCs(x, y, oriC_locations, mode_a='min', mode_b='max', window_size=window_size)
        false_order = get_false_order(sequence, x, y, oriC_locations, mode_a='min', mode_b='max', window_size=window_size)
    elif option == 2:
        oriC_locations = sort_oriCs(x, gc, oriC_locations, mode_a='min', mode_b='min', window_size=window_size)
        false_order = get_false_order(sequence, x, gc, oriC_locations, mode_a='min', mode_b='min', window_size=window_size)
    elif option == 3:
        oriC_locations = sort_oriCs(y, gc, oriC_locations, mode_a='max', mode_b='min', window_size=window_size)
        false_order = get_false_order(sequence, y, gc, oriC_locations, mode_a='max', mode_b='min', window_size=window_size)

    # Add penalties
    n_penalties = (n/len(sequence)) * len(oriC_locations)
    o_penalty   = option
    d_penalties = get_dist_penalty(len(sequence), matched_peaks)
    oriC_ranges = get_oriC_ranges(len(sequence), oriC_locations, range_size=oriC_size)

    preferred_oriC_properties = {
        'name'          : name,
        'oriC_edges'    : oriC_ranges,
        'oriC_middles'  : oriC_locations,
        'z_curve'       : (x, y, z),
        'gc_skew'       : gc,
        'n_penalty'     : n_penalties,
        'o_penalty'     : o_penalty,
        'd_penalty'     : d_penalties,
        'false_order'   : false_order,
        'seq_size'      : len(sequence),
        'gc_conc'       : ( sequence.count('G') + sequence.count('C') ) / len(sequence)
    }

    # Getting all options
    all_oriCs_list = curve_combinations( (x, y, gc), (peaks_x, peaks_y, peaks_gc), (peak_windows_x, peak_windows_y, peak_windows_gc) )

    all_oriCs_dict = {
        'xy_oriCs' : all_oriCs_list[0],
        'xgc_oriCs': all_oriCs_list[1],
        'ygc_oriCs': all_oriCs_list[2]
    }

    return preferred_oriC_properties, all_oriCs_dict


if __name__ == '__main__':
    # oriC in min of x (Purine vs. Pyrimidine)
    # oriC in max of y (Amino vs Keto)

    # # For testing a small folder
    # lengths = []
    # for fasta in os.listdir('./worst_fastas/worsened_most_v1_to_v2'):
    #     file = os.path.join('worst_fastas', 'worsened_most_v1_to_v2', fasta)
    #     # _, seq = read_FASTA(file)
    #     # lengths.append(len(seq))
    #     preferred_properties, all_oriCs = find_oriCs(file)
    
    #     name    = preferred_properties['name']
    #     Z_curve = preferred_properties['z_curve']
    #     GC_skew = preferred_properties['gc_skew']

    #     print(name)
    #     print('QoP  :', preferred_properties['n_penalty'], preferred_properties['o_penalty'], preferred_properties['d_penalty'], preferred_properties['false_order'])
    #     print('oriCs:', preferred_properties['oriC_edges'])

    #     for key in all_oriCs.keys():
    #         print(key, all_oriCs[key])

    #     pf.plot_Z_curve_2D(list(Z_curve[:2]) + [GC_skew], [preferred_properties['oriC_middles']]*3, name)
    #     # pf.plot_skew(GC_skew, [preferred_properties['oriC_middles']], name)
    #     # pf.plot_Z_curve_3D(Z_curve, name)

    # For Testing single files
    preferred_properties, all_oriCs = find_oriCs('./worst_fastas/worsened_most_v1_to_v2/NC_010175.fasta')

    name    = preferred_properties['name']
    Z_curve = preferred_properties['z_curve']
    GC_skew = preferred_properties['gc_skew']

    print(name)
    print('QoP  :', preferred_properties['n_penalty'], preferred_properties['o_penalty'], preferred_properties['d_penalty'], preferred_properties['false_order'])
    print('oriCs:', preferred_properties['oriC_edges'])

    for key in all_oriCs.keys():
        print(key, all_oriCs[key])

    pf.plot_Z_curve_2D(list(Z_curve[:2]) + [GC_skew], [preferred_properties['oriC_middles']]*3, name)
    # pf.plot_skew(GC_skew, [properties['oriC_middles']], name)
    # pf.plot_Z_curve_3D(Z_curve, name)
