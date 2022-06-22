# Libraries
import os, warnings, time
import scipy.signal as sp
import numpy as np

from itertools import combinations, product
from typing import Union, Tuple

# Self-made modules
from peak import Peak
import functions as fc
import plotting_functions as pf

# Set cwd to location of this script
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )


def calc_disparities(seq: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Z-curve and GC-skew calculation and k-mer indexing. In one function so only one iteration of the sequence is necessary.\n
    Parameters:
        `seq`        : string DNA sequence
        `k`          : length of k-mer. Should be same length as every dnaa-box
        `dnaa_boxes` : set of dnaa-box regions.\n
    Return:
        `x`, `y`, `z`, `gc` : 1D-np.arrays of the three Z-curve components and 1D-np.array of the GC-skew
        `dnaa_dict`         : Dictionary of starting indexes of dnaa-boxes in `seq`
    """
    x, y, z, gc = [], [], [], []
    a, c, t, g  = 0, 0, 0, 0

    seq_len  = len(seq)

    for i in range(seq_len):
        base = seq[i]
        if base == "A": a +=1
        elif base == "C": c +=1
        elif base == "G": g +=1
        elif base == "T": t +=1

        gc.append(g - c)
        x.append( (a + g) - (c + t) ) # Purine vs Pyrimidine
        y.append( (a + c) - (g + t) ) # Amino vs Keto
        z.append( (a + t) - (c + g) ) # Weak vs Strong Hydrogen Bonds

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
    - `curve`          : 1D-np.array
    - `peaks`          : list of Peaks of curve
    - `mode`           : 'max'|'min'. Which type of extreme do you want to find?\n
    Return:
    - `accepted_peaks` : peaks that passed both filters
    """
    rejected_peaks = []
    for peak_i, peak_j in combinations(peaks, 2):
        # Filter 1: Check if any windows intersecting the window of peak i
        if peak_i.intersecting_windows(peak_j):
            # Either peaks at the beginning and end of the DNA sequence or a peak found by sp.find_peaks() that is very close to the global min/max 
            if mode == 'max': reject_middle = np.where( curve == min(curve[peak_i.middle], curve[peak_j.middle]) )[0].tolist()
            elif mode == 'min': reject_middle = np.where( curve == max(curve[peak_i.middle], curve[peak_j.middle]) )[0].tolist()
            rejected_peaks.append(peak_i) if peak_i.middle in reject_middle else rejected_peaks.append(peak_j)

    for peak in peaks:
        # Filter 2: Check if peaks are actually the extreme in their windows
        if peak.split:
            a, b = (peak.five_side, len(curve)-1), (0, peak.three_side)
            comparator_win = a if (peak.middle >= a[0] and peak.middle <= a[1]) else b
        else:
            comparator_win = (peak.five_side, peak.three_side)

        if mode == 'max' and np.max( curve[comparator_win[0]:comparator_win[1]] ) > curve[peak.middle]:
                rejected_peaks.append(peak)
        elif mode == 'min' and np.min( curve[comparator_win[0]:comparator_win[1]] ) < curve[peak.middle]:
                rejected_peaks.append(peak)

    # Create list of peaks that passed both filters
    rejected_peaks = set(rejected_peaks)
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
    peaks_to_merge = Peak.get_peaks_to_merge(accepted_peaks)

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
        oriC_locations_list.append( [Peak(Peak.get_middle(matches[0], matches[1]), len(curves_list[0]), peaks_list[0][0].window_size) for matches in matched_peaks] )
    return oriC_locations_list


def merge_oriCs(curve_size: int, groups: list, window_size: int = 500) -> Tuple[list, list]:
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


def process_gene_loc_info(current_oriCs: list, current_occurances: list, matrix: np.ndarray) -> Tuple[list, list]:
    closest_counts = [col.argmin() for col in matrix.T]
    best_idx = max(set(closest_counts), key = closest_counts.count)

    new_oriCs = [current_oriCs[best_idx]] + current_oriCs[:best_idx] + current_oriCs[best_idx+1:]
    new_occurances = [current_occurances[best_idx]] + current_occurances[:best_idx] + current_occurances[best_idx+1:]
    # new_occurances[0] = 0.5 + new_occurances[0] if new_occurances[0] <= 0.5 else 1
    return new_oriCs, new_occurances


def find_oriCs(genome_fasta: str = None, genes_fasta: str = None, use_gene_info: bool = True, accession: str = None, email: str = None, api_key: str = None) -> dict:
    '''
    VERSION 4\n
    Locates potential oriC based on Z-curve and GC-skew analysis.
    Three window_sizes are used: 1, 3 and 5 % of the total genome length. The oriCs that were found by most combinations, get used.
    If `use_gene_info`, then the location of the DnaA and DnaN genes will be considered in the ranking of the found oriCs. The oriC closest to the genes gets promoted.

    This function either reads a given FASTA file or fetches a FASTA of a given accession directly from the NCBI database and calculates its oriC.

    Parameters:
    - `genome_fasta`  : FASTA-file with circular bacterial DNA
    - `genes_fasta`   : FASTA-file with gene info in the same format as when acquired using `E-Utils(db='nuccore', rettype='fasta_cds_na')`
    - `use_gene_info` : T|F whether to use gene_info for oriC determination. Ignored if genes_fasta is provided.
    - `accession`     : Accession number of the sequence to fetch
    - `email`         : Email adress of your NCBI account
    - `api_key`       : API Key for downloading from the NCBI database (E-Utils). Optional (as of 2022/06/08), but necessary if fetching at 10 requests per second or more.

    Return:
    - `properties` : Dictionary with oriC properties
    '''
    # NOTE: Potential further improvements:
    #   - AT% measure of oriC: should be characteristically higher than the rest of the genome.
    #   V Essential gene proximity: genes like DnaA are very likely to be close to or on the oriC
    #   - Check if DNA is circular

    # Error handling
    if genome_fasta is not None and not (genome_fasta[-5:] == 'fasta' or genome_fasta[-3:] == 'fna'):
        raise ValueError('File does not have extension \'.fasta\' or \'.fna\'. Must be FASTA-file.')
    if genome_fasta is None and genes_fasta is None and accession is None:
        raise ValueError('Did not provide files to read or accession to fetch.')
    if (accession is not None or (genes_fasta is None and use_gene_info)) and email is None:
        raise ValueError('Did not provide a email adress for fetching the accession.\n\tCreate an NCBI account at: https://www.ncbi.nlm.nih.gov/\n\tCreate an API_key at: https://www.ncbi.nlm.nih.gov/account/settings/')
    if genome_fasta is not None and accession is not None:
        warnings.warn('Provided both a fasta to read and an accession to fetch. Will ignore given accession and use accession from fasta.')

    # Sequence fetching and reading
    seq_handle = fc.fetch_file(accession, email, api_key, 'fasta') if genome_fasta is None else genome_fasta
    _accession, sequence = fc.read_FASTA(seq_handle)

    start_1 = time.time()
    # Analysing sequence properties
    dnaa_boxes = fc.get_dnaa_boxes()
    x, y, z, gc = calc_disparities(sequence)
    # Finding potential oriCs
    windows = [0.01, 0.03, 0.05]
    peaks   = []
    for fraction in windows:
        window_size = int(len(sequence) * fraction)
        peaks_x  = process_array(x , mode='min', window_size=window_size)
        peaks_y  = process_array(y , mode='max', window_size=window_size)
        peaks_gc = process_array(gc, mode='min', window_size=window_size)
        peaks.extend( [j for i in curve_combinations( (x, y, gc), (peaks_x, peaks_y, peaks_gc) ) for j in i] )

    # Finding connected components in undirected graph with a depth-first search
    matrix_pot_oriCs = fc.get_adj_mat(peaks)
    connected_groups = fc.get_connected_groups(peaks, matrix_pot_oriCs, int(len(sequence)*windows[-1]))

    # Merge potential oriCs based on their groups
    oriCs, occurances = merge_oriCs(len(sequence), connected_groups, window_size=int(len(sequence)*windows[-1]))
    checkpoint = time.time() - start_1
    if use_gene_info:
        # Gene info fetching and reading
        gene_handle = fc.fetch_file(_accession, email, api_key, 'fasta_cds_na') if genes_fasta is None else genes_fasta
        start_2 = time.time()
        genes_of_interest = ['dnaA', 'dnaN'] # 'gidA', 'parA', 'hemE' # not sure if these are proper yet
        genes_dict = fc.read_gene_info(gene_handle, genes_of_interest)
        del gene_handle

        if len(genes_dict.keys()) != 0:
            # Check oriC closest to genes of interest
            gene_locations = fc.extract_locations(len(sequence), genes_dict)
            if gene_locations is None:
                return None
            matrix_oriCs_genes = fc.get_adj_mat(oriCs, gene_locations)
            # Rearange oriCs based on gen location information
            oriCs, occurances = process_gene_loc_info(oriCs, occurances, matrix_oriCs_genes)

    # Check false order
    false_order = any( get_false_order( sequence, i[0], oriCs, mode=i[1]) for i in ((x, 'min'), (y, 'max'), (gc, 'min')) )    

    # Final dictionary
    oriC_middles = [oriC.middle for oriC in oriCs]
    finish = time.time() - start_2 + checkpoint
    preferred_oriC_properties = {
        'name'         : _accession,
        'oriC_middles' : oriC_middles,
        'occurances'   : occurances,
        'z_curve'      : (x, y, z),
        'gc_skew'      : gc,
        'false_order'  : false_order,
        'seq_size'     : len(sequence),
        'gc_conc'      : ( sequence.count('G') + sequence.count('C') ) / len(sequence),
        'time_of_prediction' : finish # Excludes time used for downloading files
    }
    print(f'Time to get oriCs: {finish:.2f} sec')
    return preferred_oriC_properties


if __name__ == '__main__':
    email = 'zoyavanmeel@gmail.com'
    api_key = '795d705fb638507c9b2295c89cc64ee88108'

    exp_refseq = [ # Accessions that have been experimentally verified.
        'NC_000964', 'NC_002947', 'NC_003272', 'NC_003869', 'NC_003888', 'NC_005090', 'NC_006461', 'NC_007633', 'NC_000913', 'NC_003047',
        'NC_007604', 'NC_000962', 'NC_002696', 'NC_002971', 'NC_005363', 'NC_008255', 'NC_009850', 'NC_010546', 'NC_010547', 'NC_011916'
    ]

    # For testing a small folder
    for fasta in exp_refseq:
        properties = find_oriCs(accession=fasta, email=email, api_key=api_key)

        name    = properties['name']
        Z_curve = properties['z_curve']
        GC_skew = properties['gc_skew']

        print(name, properties['seq_size'], 'bp \n')
        # print('QoP  :', properties['occurances'], properties['false_order'])
        # print('oriCs:', properties['oriC_middles'], '\n')


        # pf.plot_Z_curve_2D(list(Z_curve[:2]) + [GC_skew], [properties['oriC_middles']]*3, name)
        # pf.plot_skew(GC_skew, [preferred_properties['oriC_middles']], name)
        # pf.plot_Z_curve_3D(Z_curve, name)

    # For Testing single files
    # properties = find_oriCs(genome_fasta='test_fastas\Bacillus_subtilis_168.fna', email=email, api_key=api_key)
    # name    = properties['name']
    # Z_curve = properties['z_curve']
    # GC_skew = properties['gc_skew']

    # print(name)
    # print('QoP  :', properties['occurances'], properties['false_order'])
    # print('oriCs:', properties['oriC_middles'])

    # pf.plot_Z_curve_2D(list(Z_curve[:2]) + [GC_skew], [properties['oriC_middles']]*3, name)
    # # pf.plot_skew(GC_skew, [properties['oriC_middles']], name)
    # # pf.plot_Z_curve_3D(Z_curve, name)
