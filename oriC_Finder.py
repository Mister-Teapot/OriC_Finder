# Libraries
from msilib import sequence
import os
import time

import scipy.signal as sp
import numpy as np

from itertools import combinations
from itertools import product
from typing import Union
from Bio import SeqIO, Entrez

# Self-made module
import plotting_functions as pf

# Set cwd to location of this script
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

def fetch_FASTA(accession, email, api_key):
    """Downloads the nucleotide FASTA of the given accession in temporary memory"""
    Entrez.email = email
    Entrez.api_key = api_key
    return Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")


def read_FASTA(handle) -> tuple:
    """Read a FASTA file and returns the accession and sequence of only the first sequence in the file"""
    Seq_records = SeqIO.parse(handle, 'fasta')
    Seq_obj = next(Seq_records)
    return Seq_obj.id, Seq_obj.seq


class Peak():
    def __init__(self, middle: int, seq_len: int, window_size: int):
        '''
        - five_side  : 5' side of the window
        - three_side : 3' side of the window.
        '''
        self.middle = middle
        self.seq_len = seq_len
        self.window_size = window_size
        self.split = False

        five_side = self.middle - self.window_size // 2
        three_side = self.middle + self.window_size // 2

        if five_side < 0:
            self.split = True
            five_side += self.seq_len-1
        if three_side > self.seq_len-1:
            self.split = True
            three_side -= self.seq_len-1

        self.five_side = five_side
        self.three_side = three_side

    @staticmethod
    def calc_dist(curve_size: int, a: int, b: int) -> int:
        '''Calculate the distance of self to other in bp'''
        dist_1 = max(a, b) - min(a, b)
        dist_2 = min(a, b) + curve_size-1 - max(a, b)
        return min(dist_1, dist_2)

    @staticmethod
    def _get_peaks_to_merge(peaks: list) -> list:
        """
        Get the indeces that give the same value in the curve and are in eachothers window.
        Input:
        - ``peaks``          : list of Peaks of a curve
        Return:
        - ``peaks_to_merge`` : Nested list. Each sublist contains the two Peaks that have to be merged
        """
        peaks_to_merge = []
        for peak_i, peak_j in combinations(peaks, 2):
            if peak_i.contains_point(peak_j.middle):
                peaks_to_merge.append([peak_i, peak_j])
        return peaks_to_merge

    def get_middle(self, other) -> int:
        """
        Calculate the distance between Peak a and b on circular DNA.
        Returns a new Peak.middle (int) in the middle of a and b
        """
        if not isinstance(other, Peak):
            raise ValueError(f'peak_a and/or peak_b are/is not (a) Peak object(s)')
        dist_1 = max(self.middle, other.middle) - min(self.middle, other.middle)
        dist_2 = min(self.middle, other.middle) + self.seq_len-1 - max(self.middle, other.middle)
        if dist_2 < dist_1:
            merged_idx = min(self.middle, other.middle) - dist_2//2
            if merged_idx < 0:
                merged_idx = self.seq_len-1 + merged_idx
        else:
            merged_idx = min(self.middle, other.middle) + dist_1//2
        return merged_idx

    def intersecting_windows(self, other) -> bool:
        '''T|F wether self's window intersects with other's window'''
        if not isinstance(other, Peak):
            raise ValueError(f'other is a {type(other)}, must be a Peak object')
        # f = five_side, t = three_side
        f_t, t_f = Peak.calc_dist(self.seq_len, self.five_side, other.three_side), Peak.calc_dist(self.seq_len, self.three_side, other.five_side)
        f_f, t_t = Peak.calc_dist(self.seq_len, self.five_side, other.five_side), Peak.calc_dist(self.seq_len, self.three_side, other.three_side)
        a, b = f_t <= self.window_size // 2 or t_f <= self.window_size // 2, f_f <= self.window_size // 2 and t_t <= self.window_size // 2
        return a or b

    def contains_point(self, point: int) -> bool:
        '''T|F wether a point is within a Peak's window'''
        five = Peak.calc_dist(self.seq_len, self.five_side, point)
        three = Peak.calc_dist(self.seq_len, self.three_side, point)
        return (five <= self.window_size and three <= self.window_size)

    def __repr__(self):
        return f'Peak(middle={self.middle}, window_size={self.window_size})'

    def __str__(self):
        return str(self.middle)

    # Undefined __eq__ to remain hashable; only added the functions I needed
    def __lt__(self, other) -> bool:
        return self.middle < other.middle

    def __gt__(self, other) -> bool:
        return self.middle > other.middle

    def __add__(self, other) -> Union["Peak", int, float]:
        if isinstance(other, Peak):
            return Peak(self.middle + other.middle, self.seq_len, self.window_size)
        elif isinstance(other, int) or isinstance(other, float):
            return self.middle + other

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other) -> Union["Peak", int, float]:
        if isinstance(other, Peak):
            return Peak(self.middle - other.middle, self.seq_len, self.window_size)
        elif isinstance(other, int) or isinstance(other, float):
            return self.middle - other

    def __rsub__(self, other):
        return self.__sub__(other)


def calc_disparities(seq: str) -> tuple:
    """
    Z-curve, GC-skew, and N-count calculation. In one function so only one iteration of the sequence is necessary.\n
    Input:
        ``seq`` : string DNA sequence\n
    Return:
        ``x``, ``y``, ``z``, ``gc`` : 1D-np.arrays of the three Z-curve components and 1D-np.array of the GC-skew
    """
    x, y, z, gc = [], [], [], []
    a, c, t, g = 0, 0, 0, 0

    for base in seq:
        if base == "A":
            a +=1                  
        elif base == "C":
            c +=1
        elif base == "G":
            g +=1
        elif base == "T":
            t +=1 

        gc.append(g - c)
        x.append( (a + g)-(c + t) ) # Purine vs Pyrimidine
        y.append( (a + c)-(g + t) ) # Amino vs Keto
        z.append( (a + t)-(c + g) ) # Weak vs Strong Hydrogen Bonds
    return np.asarray(x), np.asarray(y), np.asarray(z), np.asarray(gc)


def detect_peaks(curve: np.ndarray) -> np.ndarray:
    """Calculates peaks of 1D-np.array and returns its indeces."""
    maxima, _ = sp.find_peaks( curve, distance=len(curve)//12)
    maxima    = np.append(maxima, curve.argmax())
    minima, _ = sp.find_peaks( np.negative(curve), distance=len(curve)//12)
    minima    = np.append(minima, curve.argmin())
    return np.unique(np.concatenate( (maxima, minima), axis=0))


def filter_peaks(curve: np.ndarray, peaks: list, mode: str = 'max') -> list:
    """
    Filters the given peaks based on the type of extreme it should be and the area around the peak in two steps.
    - Filter 1: Check if any windows intersect the window of another peak (circular DNA).
    - Filter 2: Check if peaks are actually the extreme in their windows.\n
    Input:
    - ``curve``          : 1D-np.array
    - ``peaks``          : list of Peaks of curve
    - ``mode``           : 'max'|'min'. Which type of extreme do you want to find?\n
    Return:
    - ``accepted_peaks`` : peaks that passed both filters
    """
    rejected_peaks = []
    for peak_i, peak_j in combinations(peaks, 2):
        # Filter 1: Check if any windows intersecting the window of peak i
        if peak_i.intersecting_windows(peak_j):
            # Either peaks at the beginning and end of the DNA sequence or a peak found by sp.find_peaks() that is very close to the global min/max 
            if mode == 'max':
                reject_middle = np.where( curve == min(curve[peak_i.middle], curve[peak_j.middle]) )[0].tolist()
            elif mode == 'min':
                reject_middle = np.where( curve == max(curve[peak_i.middle], curve[peak_j.middle]) )[0].tolist()
            rejected_peaks.append(peak_i) if peak_i.middle in reject_middle else rejected_peaks.append(peak_j)

    for peak in peaks:
        # Filter 2: Check if peaks are actually the extreme in their windows
        if peak.split:
            a = (peak.five_side, len(curve)-1)
            b = (0, peak.three_side)
            comparator_win = a if (peak.middle >= a[0] and peak.middle <= a[1]) else b
        else:
            comparator_win = (peak.five_side, peak.three_side)

        if mode == 'max' and np.max( curve[comparator_win[0]:comparator_win[1]] ) > curve[peak.middle]:
                rejected_peaks.append(peak)
        elif mode == 'min' and np.min( curve[comparator_win[0]:comparator_win[1]] ) < curve[peak.middle]:
                rejected_peaks.append(peak)

    # Create list of peaks that passed both filters
    rejected_peaks = list(set(rejected_peaks))
    accepted_peaks = [x for x in peaks if x not in rejected_peaks]
    return accepted_peaks


def match_peaks(peaks_x: list, peaks_y: list) -> list:
    """Nested list of peaks from x that line up with peaks from y."""
    matched_peaks = []
    for peak_x, peak_y in product(peaks_x, peaks_y):
        if peak_x.intersecting_windows(peak_y):
            matched_peaks.append( (peak_x, peak_y) )
    return matched_peaks


def process_array(curve: np.ndarray, mode: str = 'max', window_size: int = 500) -> list:
    """
    Runs the given 1D-array (curve) through all processing functions for oriC identification.
    Returns its peaks and windows.
    """
    init_peaks = [Peak(peak, len(curve), window_size) for peak in detect_peaks(curve)]
    accepted_peaks = filter_peaks(curve, init_peaks, mode=mode)
    peaks_to_merge = Peak._get_peaks_to_merge(accepted_peaks)

    single_peaks = [x for x in accepted_peaks if not any(x in y for y in peaks_to_merge)]
    merged_peaks = [Peak(to_merge[0].get_middle(to_merge[1]),len(curve), window_size) for to_merge in peaks_to_merge]
    return single_peaks + merged_peaks


def get_false_order(seq: str, curve: np.ndarray, oriC_locations: list, mode: str = 'max') -> bool:
    '''T|F, whether the order of the oriCs could be wrong due to 'N' bases in area around oriC'''
    false_orders = []
    if len(oriC_locations) > 1:
        n_per_oriC = []
        for oriC in oriC_locations:
            if oriC.split:
                n = seq[:oriC.three_side].count('N') + seq[oriC.five_side:].count('N')
            else:
                n = seq[oriC.three_side:oriC.five_side].count('N')
            n_per_oriC.append(n)

        # NOTE: only checks the oriC in top position against the other oriCs, not every combination
        # This looks at cases where there are two peaks of similar heights/depths that are quite far apart.
        # e.g. if the x-curve was W-shaped, making it so you have two very similar minima. (or y-curve like M)
        for i in range(1, len(n_per_oriC)):
            if mode == 'max':
                false_orders.append( curve[oriC_locations[0].middle] - n_per_oriC[0] <= curve[oriC_locations[i].middle] + n_per_oriC[i] )
            elif mode == 'min':
                false_orders.append( curve[oriC_locations[0].middle] + n_per_oriC[0] >= curve[oriC_locations[i].middle] - n_per_oriC[i] )
    return any(false_orders)


def curve_combinations(curves_list: Union[list, tuple], peaks_list: Union[list, tuple]) -> list:
    '''Get every matched_peaks combination for x, y, and gc.'''
    oriC_locations_list = []
    for peaks_i, peaks_j in combinations(peaks_list, 2):
        matched_peaks  = match_peaks(peaks_i, peaks_j)
        oriC_locations_list.append( [Peak(matches[0].get_middle(matches[1]), len(curves_list[0]), peaks_list[0][0].window_size) for matches in matched_peaks] )
    return oriC_locations_list


def get_adj_mat(peaks: list) -> np.ndarray:
    '''Gets adjacency matrix for given peaks'''
    adj_mat = np.zeros((len(peaks), len(peaks)))
    for (i_a, a), (i_b, b) in combinations(enumerate(peaks), r=2):
        dist = Peak.calc_dist(a.seq_len, a.middle, b.middle)
        adj_mat[i_a, i_b] = dist
        adj_mat[i_b, i_a] = dist
    return adj_mat


def get_connected_groups(peaks: list, adj_mat: np.ndarray, threshold: int) -> list:
    '''Recursively find connected groups in an undirected graph'''
    visited = [False] * len(peaks)
    connected_groups_idx = []
    for i in range(len(peaks)):
        if not visited[i]:
            group = []
            _, _, visited, group, _ = __DFS_recurse(i, adj_mat, visited, group, threshold=threshold)
            connected_groups_idx.append(group)
    connected_groups_vals = [ [peaks[i] for i in idx_group] for idx_group in connected_groups_idx ]
    return connected_groups_vals


def __DFS_recurse(idx, adj_mat, visited, connected_list, threshold):
    visited[idx] = True
    connected_list.append(idx)
    for i in range(len(visited)):
        if i == idx:
            continue
        elif adj_mat[i][idx] <= threshold and not visited[i]:
            _, _, visited, connected_list, _ = __DFS_recurse(i, adj_mat,visited, connected_list, threshold)
    return idx, adj_mat, visited, connected_list, threshold


def merge_oriCs(curve_size: int, groups: list, window_size: int = 500) -> tuple:
    '''Finds the average index of a group and returns those values. groups is a nested-list'''
    mutable = sorted( groups, key=lambda x:len(x), reverse=True )
    oriCs, occurances = [], []
    for group in mutable:
        group.sort()
        for i in range(len(group)):
            if (group[-1] - group[i]).middle >= (curve_size-1)/2:
                group[i].middle += curve_size-1
        avg_val = sum(group)//len(group)
        if avg_val > curve_size-1:
            avg_val -= curve_size-1
        oriCs.append(Peak(avg_val, curve_size, window_size))
        occurances.append(len(group))

    total_pot_oriCs = len( [y for x in mutable for y in x] )
    occurances = [x/total_pot_oriCs for x in occurances]
    return oriCs, occurances


def find_oriCs(filename: str = None, accession: str = None, email: str = None, API_key: str = None) -> dict:
    '''
    VERSION 4\n
    Locates potential oriC based on Z-curve and GC-skew analysis.
    Three window_sizes are used: 1, 3 and 5 % of the total genome length. The oriCs that were found by most combinations, get returned.
    TODO: Add DnaA-box analysis

    This function either reads a given FASTA file or fetches a FASTA of a given accession directly from the NCBI database.

    Parameters:
    - ``filename`` : FASTA file with circular bacterial DNA
    - ``accession``: Accession number of the sequence to fetch
    - ``email``    : Email adress of your NCBI account
    - ``API_key``  : API Key for downloading from the NCBI database (E-Utils).

    Return:
    - ``properties``  : Dictionary with oriC properties
    '''
    # NOTE: Potential further improvements:
    #   - AT% measure of oriC: should be characteristically higher than the rest of the genome.
    #   - Essential gene proximity: genes like DnaA are very likely to be close to or on the oriC

    # Error handling
    if filename is None and accession is None:
        raise ValueError('Did not provide a FASTA to read or accession to fetch.')
    if filename is not None and accession is not None:
        raise ValueError('Provided both a file to read and an accession to fetch, choose one.')
    if accession is not None and (email is None or API_key is None):
        raise ValueError('Did not provide a email adress or API key for fetching the sequence.\n\tCreate an NCBI account at: https://www.ncbi.nlm.nih.gov/\n\tCreate an API_key at: https://www.ncbi.nlm.nih.gov/account/settings/')

    # Fetching and reading
    handle = fetch_FASTA(accession, email, API_key) if accession is not None else filename
    accession, seq = read_FASTA(handle)
    del handle # either the path to a file or an entire FASTA file

    # Analysing sequence properties
    x, y, z, gc = calc_disparities(seq)

    windows = [0.01, 0.03, 0.05]
    peaks   = []
    for fraction in windows:
        window_size = int(len(seq) * fraction)
        peaks_x  = process_array(x , mode='min', window_size=window_size)
        peaks_y  = process_array(y , mode='max', window_size=window_size)
        peaks_gc = process_array(gc, mode='min', window_size=window_size)
        peaks.extend( [y for x in curve_combinations( (x, y, gc), (peaks_x, peaks_y, peaks_gc) ) for y in x] )

    # Connected components in undirected graph problem
    matrix = get_adj_mat(peaks)

    # Depth-First Search
    connected_groups = get_connected_groups(peaks, matrix, int(len(seq)*windows[-1]))

    # Remove potential oriC if it was not matched to any other.
    connected_groups = [x for x in connected_groups if len(x) > 1]

    # Get oriCs
    oriCs, occurances = merge_oriCs(len(seq), connected_groups, window_size=int(len(seq)*windows[-1]))

    # Check false order
    false_order = any( get_false_order( seq, i[0], oriCs, mode=i[1]) for i in ((x, 'min'), (y, 'max'), (gc, 'min')) )        

    # Final dictionary
    preferred_oriC_properties = {
        'name'         : accession,
        'oriC_middles' : [oriC.middle for oriC in oriCs],
        'occurances'   : occurances,
        'z_curve'      : (x, y, z),
        'gc_skew'      : gc,
        'false_order'  : false_order,
        'seq_size'     : len(seq),
        'gc_conc'      : ( seq.count('G') + seq.count('C') ) / len(seq)
    }

    return preferred_oriC_properties


if __name__ == '__main__':
    email = 'zoyavanmeel@gmail.com'
    API_key='795d705fb638507c9b2295c89cc64ee88108'
    # For testing a small folder
    # for fasta in os.listdir('./test_fastas'):
    #     file = os.path.join('test_fastas', fasta)
    #     properties = find_oriCs(file)

    #     name    = properties['name']
    #     Z_curve = properties['z_curve']
    #     GC_skew = properties['gc_skew']

    #     print(name)
    #     print('QoP  :', properties['occurances'], properties['false_order'])
    #     print('oriCs:', properties['oriC_middles'])


    #     pf.plot_Z_curve_2D(list(Z_curve[:2]) + [GC_skew], [properties['oriC_middles']]*3, name)
    #     # pf.plot_skew(GC_skew, [preferred_properties['oriC_middles']], name)
    #     # pf.plot_Z_curve_3D(Z_curve, name)

    # For Testing single files
    properties = find_oriCs(filename='./test_fastas/Bacillus_subtilis_168.fna')
    name    = properties['name']
    Z_curve = properties['z_curve']
    GC_skew = properties['gc_skew']

    print(name)
    print('QoP  :', properties['occurances'], properties['false_order'])
    print('oriCs:', properties['oriC_middles'])

    pf.plot_Z_curve_2D(list(Z_curve[:2]) + [GC_skew], [properties['oriC_middles']]*3, name)
    # pf.plot_skew(GC_skew, [properties['oriC_middles']], name)
    # pf.plot_Z_curve_3D(Z_curve, name)
