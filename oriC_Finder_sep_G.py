# Libraries
import os, warnings, time
import scipy.signal as sp
import numpy as np
import joblib

from sklearn import svm
from itertools import combinations, product
from typing import Union, Tuple

# Self-made modules
from peak import Peak
import functions as fc
import plotting_functions as pf

# Set cwd to location of this script
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )


def calc_disparities(seq: str, k: int, dnaa_boxes: set) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    '''
    Z-curve and GC-skew calculation and k-mer indexing. In one function so only one iteration of the sequence is necessary.\n
    Parameters:
        `seq`        : string DNA sequence
        `k`          : length of k-mer. Should be same length as every dnaa-box
        `dnaa_boxes` : set of dnaa-box regions.\n
    Return:
        `x`, `y`, `z`, `gc` : 1D-np.arrays of the three Z-curve components and 1D-np.array of the GC-skew
        `dnaa_dict`         : Dictionary of starting indexes of dnaa-boxes in `seq`
    '''
    x, y, z, gc = [], [], [], []
    a, c, t, g  = 0, 0, 0, 0

    raw_dict = {}
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

        kmer = seq[i:i+k] if i <= seq_len - k else seq[i:] + seq[:k-(seq_len-i)]
        mid = i+5 if i+5 <= seq_len else i+5 - seq_len
        try: raw_dict[kmer].append(mid)
        except KeyError: raw_dict[kmer] = [mid]

    keys = dnaa_boxes & raw_dict.keys() # & assumes .keys() as a set (which it should as dict keys are unique). .intersection() assumes .keys as a list and set.intersection(list) has a worse time-complexity. https://wiki.python.org/moin/TimeComplexity
    dnaa_dict = {key : raw_dict[key] for key in keys}
    del raw_dict
    return np.asarray(x), np.asarray(y), np.asarray(z), np.asarray(gc), dnaa_dict


def detect_peaks(curve: np.ndarray) -> np.ndarray:
    '''Calculates peaks of 1D-np.array and returns its indeces.'''
    maxima, _ = sp.find_peaks( curve, distance=len(curve)//12)
    maxima    = np.append(maxima, curve.argmax())
    minima, _ = sp.find_peaks( np.negative(curve), distance=len(curve)//12)
    minima    = np.append(minima, curve.argmin())
    return np.unique(np.concatenate( (maxima, minima), axis=0))


def filter_peaks(curve: np.ndarray, peaks: list, mode: str) -> list:
    '''
    Filters the given peaks based on the type of extreme it should be and the area around the peak in two steps.
    - Filter 1: Check if any windows intersect the window of another peak (circular DNA).
    - Filter 2: Check if peaks are actually the extreme in their windows.\n
    Input:
    - `curve`          : 1D-np.array
    - `peaks`          : list of Peaks of curve
    - `mode`           : 'max'|'min'. Which type of extreme do you want to find?\n
    Return:
    - `accepted_peaks` : peaks that passed both filters
    '''
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
    '''Nested list of peaks from x that line up with peaks from y.'''
    matched_peaks = []
    for peak_x, peak_y in product(peaks_x, peaks_y):
        if peak_x.intersecting_windows(peak_y):
            matched_peaks.append( (peak_x, peak_y) )
    return matched_peaks


def get_false_order(seq: str, curve: np.ndarray, oriC_locations: list, mode: str) -> bool:
    '''T|F, whether the order of the oriCs could be wrong due to 'N' bases in area around oriC. mode: min|max'''
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


def merge_oriCs(curve_size: int, groups: list, window_size: int) -> Tuple[list, list]:
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


def process_array(curve: np.ndarray, mode: str, window_size: int) -> list:
    '''Runs the given 1D-array (curve) through all processing functions for oriC identification. Returns its peaks'''
    init_peaks = [Peak(peak, len(curve), window_size) for peak in detect_peaks(curve)]
    accepted_peaks = filter_peaks(curve, init_peaks, mode=mode)
    peaks_to_merge = Peak.get_peaks_to_merge(accepted_peaks)

    single_peaks = [x for x in accepted_peaks if not any(x in y for y in peaks_to_merge)]
    merged_peaks = [Peak(to_merge[0].get_middle(to_merge[1]),len(curve), window_size) for to_merge in peaks_to_merge]
    return single_peaks + merged_peaks


def get_occurances_gene_loc_info(matrix: np.ndarray) -> list:
    '''Process the location of the genes of interest and rank the potential oriCs based on how close they are to these genes'''
    if np.min(matrix) == np.max(matrix):
        return [[1]] * matrix.shape[1]
    All_G_occurances = []
    for col in range(matrix.shape[1]):
        if np.max(matrix[:,col]) == np.min(matrix[:,col]):
            G = np.asarray([1])
        else:
            G = (matrix[:,col] - np.min(matrix[:,col])) / (np.max(matrix[:,col]) - np.min(matrix[:,col]))
        All_G_occurances.append(G.tolist())
    return All_G_occurances


def get_occurances_box_loc_info(current_oriCs: list, kmer_dict: dict) -> list:
    '''Process the location of dnaa_boxes and rank potential oriCs based on most surrounding dnaa_boxes.'''
    contains_boxes = []
    all_pos = [pos for pos_list in kmer_dict.values() for pos in pos_list]
    for oriC in current_oriCs:
        count = 0
        for pos in all_pos:
            if oriC.contains_point(pos):
                count += 1
        contains_boxes.append(count)
    return [x/sum(contains_boxes) if sum(contains_boxes) != 0 else 0 for x in contains_boxes]


def find_oriCs(
    genome_fasta: str = None,
    genes_fasta: str = None,
    custom_dnaa_boxes: list = None,
    windows: list = None,
    threshold: float = None,
    max_mismatches: int = 0,
    model: svm.SVC = None,
    accession: str = None,
    email: str = None,
    api_key: str = None
) -> dict:
    '''
    VERSION 6\n
    Locates potential oriC on circular bacterial chromosomes based on Z-curve and GC-skew analysis and dnaa box analysis.
    Three window_sizes are used: 1, 3 and 5 % of the total genome length. The oriCs that were found by most combinations, get used.
    The location of the DnaA and DnaN genes will be considered in the ranking of the found oriCs. The oriC closest to the genes gets promoted.

    This function either reads a given FASTA and genes_fasta or fetches them using a given accession directly from the NCBI database and calculates its oriC.

    Parameters:
    - `genome_fasta`        : FASTA-file with circular bacterial DNA
    - `genes_fasta`         : FASTA-file with gene info in the same format as when acquired using `E-Utils(db='nuccore', rettype='fasta_cds_na')`
    - `custom_dnaa_boxes`   : If None, will use the [consensus DnaA-box](https://doi.org/10.1093/bib/bbn031): `TTAT(A|T|C|G)CACA`.
                              Else, provide a list of 9 base strings. Regex **NOT** supported for all-unique box generation.
                              Example input: `['AAAAAAAAA', 'TTTTTTTTT']`.
    - `max_mismatches`      : Maximum allowed mismatches before a 9-mer is considered to fit the dnaa_box. Recommended: 0; recommended max: 2.
    - `accession`           : Accession number of the sequence to fetch
    - `email`               : Email adress of your NCBI account
    - `api_key`             : API Key for downloading from the NCBI database (E-Utils). Optional (as of 2022/06/22), but necessary if fetching at 10 requests per second or more.

    Return:
    - `properties` : Dictionary with oriC properties
    '''
    # NOTE: Potential further improvements:
    #   - AT% measure of oriC: should be characteristically higher than the rest of the genome.

    # Clarification:
    #   Z_oriCs: order of found oriCs according to only the Z-curve
    #   G_oriCs: order of Z_oriCs based on proximity to genes
    #   D_oriCs: order of Z_oriCs based on number of DnaA-boxes in its 5%-window.

    # Some error handling
    if genome_fasta is not None and not (genome_fasta[-5:] == 'fasta' or genome_fasta[-3:] == 'fna'):
        raise ValueError('\'genome_fasta\' does not have extension \'.fasta\' or \'.fna\'. Must be FASTA-file.')
    if genes_fasta is not None and not (genes_fasta[-5:] == 'fasta' or genes_fasta[-3:] == 'fna'):
        raise ValueError('\'genes_fasta\' does not have extension \'.fasta\' or \'.fna\'. Must be FASTA-file.')
    if genome_fasta is None and genes_fasta is None and accession is None:
        raise ValueError('Did not provide files to read or accession to fetch.')
    if genome_fasta is not None and accession is not None:
        warnings.warn('Provided both a fasta to read and an accession to fetch. Will ignore given accession and use accession from \'genome_fasta\'.')
    if accession is not None and email is None:
        raise ValueError('Did not provide a email adress for fetching the accession.\n\tCreate an NCBI account at: https://www.ncbi.nlm.nih.gov/\n\tCreate an API_key at: https://www.ncbi.nlm.nih.gov/account/settings/')
    if genome_fasta is not None and genes_fasta is None and email is None:
        raise ValueError('Only provided \'genome_fasta\'. Will have to fetch \'genes_fasta\', but you provided no \'email\'.')

    # Sequence fetching and reading
    seq_handle = fc.fetch_file(accession, email, api_key, 'fasta') if genome_fasta is None else genome_fasta
    start_1 = time.time()
    _accession, sequence = fc.read_FASTA(seq_handle)
    seq_len = len(sequence)
    gc_conc = ( sequence.count('G') + sequence.count('C') ) / seq_len

    # User-defined variables
    box_list = ['TTATACACA', 'TTATTCACA', 'TTATCCACA', 'TTATGCACA'] if custom_dnaa_boxes is None else list(custom_dnaa_boxes)
    _windows = [0.01, 0.03, 0.05] if windows is None else windows
    _threshold = int(seq_len*_windows[-1]) if threshold is None else threshold

    # Analysing sequence properties
    dnaa_boxes = fc.get_dnaa_boxes(box_list=box_list, max_mismatches=int(max_mismatches))
    start_calc = time.time()
    x, y, z, gc, kmers = calc_disparities(sequence, 9, dnaa_boxes)
    calc_disp_time = time.time()-start_calc
    del sequence

    # Finding potential oriCs based on Z-curve
    peaks   = []
    for fraction in _windows:
        window_size = int(seq_len * fraction)
        peaks_x  = process_array(x , mode='min', window_size=window_size)
        peaks_y  = process_array(y , mode='max', window_size=window_size)
        peaks_gc = process_array(gc, mode='min', window_size=window_size)
        peaks.extend( [j for i in curve_combinations( (x, y, gc), (peaks_x, peaks_y, peaks_gc) ) for j in i] )

    ## Finding connected components in undirected graph with a depth-first search to merge Z-curve oriCs
    # pf.plot_Z_curve_2D([x, y, gc], [[peak.middle for peak in peaks]]*3, ['$x_n$', '$y_n$', '$g_n$'])
    matrix_pot_oriCs      = fc.get_adj_mat(peaks)
    connected_groups      = fc.get_connected_groups(peaks, matrix_pot_oriCs, _threshold)
    Z_oriCs, Z_occurances = merge_oriCs(seq_len, connected_groups, window_size=int(seq_len*_windows[-1]))
    rank_dict = {oriC : [occurance] for oriC, occurance in zip(Z_oriCs, Z_occurances)}

    # Gene info fetching and reading
    checkpoint = time.time() - start_1
    gene_handle = fc.fetch_file(_accession, email, api_key, 'fasta_cds_na') if genes_fasta is None else genes_fasta
    start_2 = time.time()
    genes_of_interest = ['dnaA', 'dnaN'] # 'gidA', 'parA', 'hemE' # not sure if these are proper yet
    genes_dict, num_of_genes = fc.read_gene_info(gene_handle, genes_of_interest)
    read_genes_time = time.time()-start_2
    del gene_handle

    # Gene-location analysis
    if len(genes_dict.keys()) != 0:
        gene_locations = fc.extract_locations(seq_len, genes_dict)
        if gene_locations is None: return None
        matrix_oriCs_genes = fc.get_adj_mat(Z_oriCs, gene_locations)
        G_dnaA_occurances, G_dnaN_occurances = get_occurances_gene_loc_info(matrix_oriCs_genes)
    else:
        warnings.warn(f'\n\n\tAccession: {_accession}. None of the genes of interest were found in the \'genes_fasta\': {genes_of_interest}\n\tWill not use gene locations in prediction.\n')
        G_dnaA_occurances = [0] * len(Z_oriCs)
        G_dnaN_occurances = [0] * len(Z_oriCs)

    # DnaA-box analysis
    if len(kmers.keys()) != 0:
        D_occurances = get_occurances_box_loc_info(Z_oriCs, kmers)
    else:
        warnings.warn(f'\n\n\tAccession: {_accession}. No DnaA-boxes were found: {dnaa_boxes}\n\tWill not use DnaA-boxes in prediction.\n')
        D_occurances = [0] * len(Z_oriCs)

    Avg_occurances = []
    for i, key in enumerate(rank_dict.keys()):
        rank_dict[key].append(G_dnaA_occurances[i])
        rank_dict[key].append(G_dnaN_occurances[i])
        rank_dict[key].append(D_occurances[i])
        avg = sum(rank_dict[key])/len(rank_dict[key])
        rank_dict[key].append(avg)
        Avg_occurances.append(avg)

    '''
    rank_dict = {
              # Z_score     G_score     D_score   Avg_score
    'oriC_1' : [   0.7,        0.4,        0.4,       0.50],
    'oriC_2' : [   0.2,        0.5,        0.4,       0.37],
    'oriC_3' : [   0.1,        0.1,        0.2,       0.13]
    }
    '''

    # Final dictionary: sort based on Avg_occurances
    oriCs = [oriC.middle for oriC in sorted( rank_dict.keys(), key = lambda x:rank_dict[x][-1], reverse=True )]
    Z_occurances = sorted(Z_occurances, key = lambda x:Avg_occurances[Z_occurances.index(x)], reverse=True )
    G_dnaA_occurances = sorted(G_dnaA_occurances, key = lambda x:Avg_occurances[G_dnaA_occurances.index(x)], reverse=True )
    G_dnaN_occurances = sorted(G_dnaN_occurances, key = lambda x:Avg_occurances[G_dnaN_occurances.index(x)], reverse=True )
    D_occurances = sorted(D_occurances, key = lambda x:Avg_occurances[D_occurances.index(x)], reverse=True )
    Avg_occurances = sorted( Avg_occurances, reverse=True )
    total_time = time.time() - start_2 + checkpoint

    decision = None
    if model is not None:
        decision = model.decision_function(np.asarray([Z_occurances,G_dnaA_occurances, G_dnaN_occurances, D_occurances]).T,).tolist()

    oriC_properties = {
        'name'         : _accession,
        'oriC_middles' : oriCs,
        'occurances'   : Avg_occurances,
        'Z_occurances' : Z_occurances,
        'G_A_occurances' : G_dnaA_occurances,
        'G_N_occurances' : G_dnaN_occurances,
        'D_occurances' : D_occurances,
        'Prediction'   : decision,
        'z_curve'      : (x, y, z),
        'gc_skew'      : gc,
        'seq_size'     : seq_len,
        'gc_conc'      : gc_conc,
        'time_of_prediction' : (calc_disp_time, read_genes_time, total_time), # Excludes time used for downloading files
        'num_of_genes'       : num_of_genes
    }

    if 0: print(f'Time to get oriCs: {total_time:.2f} sec')
    return oriC_properties


if __name__ == '__main__':
    email = 'zoyavanmeel@gmail.com'
    api_key = '795d705fb638507c9b2295c89cc64ee88108'
    model = joblib.load('model.pkl')

    exp_refseq = [ # Accessions that have been experimentally verified.
        'NC_000964', 'NC_002947', 'NC_003272', 'NC_003869', 'NC_005090', 'NC_006461', 'NC_007633', 'NC_000913', 'NC_003047',
        'NC_007604', 'NC_000962', 'NC_002696', 'NC_002971', 'NC_005363', 'NC_008255', 'NC_009850', 'NC_010546', 'NC_010547', 'NC_011916'
    ]

    worst_exp_refseq = ['NC_007604', 'NC_010546', 'NC_002971', 'NC_000117']

    # For testing a small folder
    # for acc in exp_refseq:
    #     properties = find_oriCs(accession=acc, email=email, api_key=api_key)

    #     name    = properties['name']
    #     Z_curve = properties['z_curve']
    #     GC_skew = properties['gc_skew']

    #     print(name, properties['oriC_middles'][0])
    #     # print(name, properties['seq_size'], 'bp \n')
    #     # print('QoP  :', properties['occurances'], properties['false_order'])
    #     # print('oriCs:', properties['oriC_middles'], '\n')

    #     pf.plot_Z_curve_2D(list(Z_curve[:2]) + [GC_skew], [properties['oriC_middles']]*3, name)
    #     # pf.plot_skew(GC_skew, [preferred_properties['oriC_middles']], name)
    #     # pf.plot_Z_curve_3D(Z_curve, name)

    # For Testing single files

    # NC_016609: good example of 'harder' sequence. Can't just look for global extremes
    properties = find_oriCs(accession='NC_021985', email=email, api_key=api_key, model=None)
    name    = properties['name']
    Z_curve = properties['z_curve']
    GC_skew = properties['gc_skew']

    print(name)
    print('Len  :', properties['seq_size'])
    print('QoP  :', properties['occurances'])
    print('Pred :', properties['Prediction'])
    print('oriCs:', properties['oriC_middles'])

    pf.plot_Z_curve_2D(list(Z_curve[:2]) + [GC_skew], [properties['oriC_middles']]*3, ['$x_n$', '$y_n$', '$g_n$'])
    # pf.plot_skew(GC_skew, [properties['oriC_middles']], name)
    # pf.plot_Z_curve_3D(Z_curve, name)