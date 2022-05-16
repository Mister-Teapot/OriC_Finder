import os, sys
from ast import literal_eval
import pandas as pd
import warnings

# Pandas printing options
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 15)

# Self-made module
#   Local path
sys.path.append('../OriC_Finder')
#   Cluster path
# sys.path.append('/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/OriC_Finder/')
import oriC_Finder as of
import plotting_functions as pf

MAX_DIST = 3/100


def make_comparator_csv(Z_curve_csv, DoriC_csv, comparator_csv):
    '''Merges Z_curve_csv and DoriC_csv into one csv and/or dataframe for comparison with compare_dbs'''

    # Read CSVs
    Z_curve_df = pd.read_csv( Z_curve_csv )
    DoriC_df   = pd.read_csv( DoriC_csv )

    # Drop empty columns
    Z_curve_df.dropna(axis=1, how='all', inplace=True)

    # Drop duplicate rows
    Z_curve_df.drop_duplicates(subset=['RefSeq'], ignore_index=True, inplace=True)

    # Remove version numbers from accessions
    Z_curve_df['RefSeq'] = Z_curve_df['RefSeq'].str.extract(r'([^.]*)')

    # Merge into single dataframe
    col_names_Z = [x for x in Z_curve_df.columns if x not in ['RefSeq', 'Organism']]
    in_both_df = {x: [] for x in DoriC_df.columns.to_list() + col_names_Z}
    for i, sample in Z_curve_df.iterrows():
        if sample['RefSeq'] in DoriC_df['RefSeq'].to_list():
            for j in col_names_Z:
                in_both_df[j].append( sample[j] )
            for j in DoriC_df.columns.to_list():
                in_both_df[j].append( DoriC_df.loc[DoriC_df['RefSeq'] == sample['RefSeq'], j].values[0] ) # Series.values -> ['element'], very annoying
    
    comparator_df = pd.DataFrame(in_both_df)
    comparator_df.to_csv(comparator_csv, index=False)
    return comparator_df


def compare_dbs(df=None, csv=None, info_file_path=None, print_info=False, max_dist=30000, n=None, d=None, alt_oriC=None):
    '''
    Compares results. Results can be loaded from a pandas.DataFrame OR csv.
    info       : If str, writes a text_file to the given location with some info on the analysis.
    print_info : If True, prints the same info that was written to 'info' with print-statements for quick reading.
    max_dist   : Maximum distance two oriCs can have from each other if they want to be matched. If <= 1, will be
                 seen as a fraction of the total sequence length. If > 1, it will be seen as # of basepairs.
    '''

    # Errors and warnings
    if df is None and csv is None:
        raise KeyError('Must load a pandas.DataFrame or CSV file.')
    if df is not None and csv is not None:
        raise KeyError('Load either a pandas.DataFrame or CSV, not both.')
    if csv is not None:
        df = pd.read_csv(csv)
    if max_dist > 30000:
        warnings.warn("Warning: Raising the max_dist above the minimum distance of oriC found with the predictor. Potential overlap in DoriC to Z_oriC matching")

    # Compare predictions (Assuming DoriC as ground truth for now)
    num_of_predictions = [0, 0] # [more oriC in DoriC, more oriC by mine]
    multi_matching = []

    all_distances = {
        'RefSeq'  : [],
        'D_idx'   : [],
        'Z_idx'   : [],
        'Distance': []
    }

    total_D_oriC = 0
    total_Z_oriC = 0
    total_D_by_Z = 0

    for i, sample_original in df.iterrows():
        sample = sample_original.copy()
        sample.dropna(inplace=True)

        seq_len = sample['Sequence_length']
        D_oriC_cols = sorted( [i for i in sample.axes[0] if 'DoriC_oriC' in i], key=lambda x:x[-1] )
        if alt_oriC is None:
            Z_oriC_cols = sorted( [i for i in sample.axes[0] if 'oriC_middle' in i], key=lambda x:x[-1] )[:2] # only use top 2 oriC if there is more than one
        else:
            Z_oriC_cols = sorted( [i for i in sample.axes[0] if alt_oriC in i], key=lambda x:x[-1] )[:2]
    
        total_D_oriC += len(D_oriC_cols)
        total_Z_oriC += len(Z_oriC_cols)
    
        if len(D_oriC_cols) > len(Z_oriC_cols):
            num_of_predictions[0] += 1
        elif len(D_oriC_cols) < len(Z_oriC_cols):
            num_of_predictions[1] += 1

        # We match each Z_oriC to the closest DoriC_oriC with a max distance of 30k basepairs (like the peak_windows, but can be adjusted later).
        # If multiple Z_oriC match with the same DoriC_oriC : handle later + check vice versa too.
        D_oriC_middles = [ of.merge_peaks( seq_len, int( literal_eval(sample[D_oriC_col])[0] ), int( literal_eval(sample[D_oriC_col])[1] ) ) for D_oriC_col in D_oriC_cols ]
        Z_oriC_middles = [ int(sample[Z_oriC_col]) for Z_oriC_col in Z_oriC_cols ]

        ### Accuracy per point: how many points were properly matched
        window_size = max_dist if max_dist > 1 else int(seq_len * max_dist)
        D_oriC_windows = of.get_peak_windows(seq_len, D_oriC_middles, window_size=window_size)
        Z_oriC_windows = of.get_peak_windows(seq_len, Z_oriC_middles, window_size=window_size)

        # Match all Z_oriC to D_oriC
        all_matches = of.match_peaks(D_oriC_middles, Z_oriC_middles, D_oriC_windows, Z_oriC_windows)

        # Match only the top len(D_oriC) Z_oriC to D_oriC
        top_matches = of.match_peaks(D_oriC_middles, Z_oriC_middles[:len(D_oriC_middles)], D_oriC_windows, Z_oriC_windows[:len(D_oriC_middles)])

        # Only work with top_matches for now (it does not differ from all_matches right now)
        D_appearances = [match[0] for match in top_matches]
        Z_appearances = [match[1] for match in top_matches]
        total_D_by_Z += len(D_appearances)
        multi_matching.append( ( len(D_appearances) != len(set(D_appearances)), len(Z_appearances) != len(set(Z_appearances)) ) )

        for i, D_oriC in enumerate(D_oriC_middles):
            distances = []
            for j, Z_oriC in enumerate(Z_oriC_middles):
                if d is None and n is None:
                    dist_1 = max(D_oriC, Z_oriC) - min(D_oriC, Z_oriC)
                    dist_2 = min(D_oriC, Z_oriC) + seq_len-1 - max(D_oriC, Z_oriC)
                    distances.append( (j, min(dist_1, dist_2)) )
                elif d is not None:
                    if sample['D_penalty'][j] <= d:
                        dist_1 = max(D_oriC, Z_oriC) - min(D_oriC, Z_oriC)
                        dist_2 = min(D_oriC, Z_oriC) + seq_len-1 - max(D_oriC, Z_oriC)
                        distances.append( (j, min(dist_1, dist_2)) )
                elif n is not None:
                    if sample['N_penalty'] <= n:
                        dist_1 = max(D_oriC, Z_oriC) - min(D_oriC, Z_oriC)
                        dist_2 = min(D_oriC, Z_oriC) + seq_len-1 - max(D_oriC, Z_oriC)
                        distances.append( (j, min(dist_1, dist_2)) )
            if len(distances) != 0:
                sorted_dists = min(distances, key=lambda x:abs(x[1]))

                all_distances['RefSeq'].append(sample['RefSeq'])
                all_distances['D_idx'].append(i)
                all_distances['Z_idx'].append(sorted_dists[0])
                all_distances['Distance'].append(sorted_dists[1])

        ### Accuracy per sample: how many samples were properly matched:

    # Information
    info_lines = ['INFORMATION']
    info_lines.append( f'There were {df.shape[0]} organisms.' )
    info_lines.append( f'\t{df.shape[0] - sum(num_of_predictions)} times both databases predicted the same amount of oriCs.' )
    info_lines.append( f'\t{num_of_predictions[0]} times DoriC predicted more oriCs.' )
    info_lines.append( f'\t{num_of_predictions[1]} times mine predicted more oriCs.\n' )

    did_not = 'did not' if all_matches == top_matches else ''
    info_lines.append( f'If a DoriC sample had 2 oriC and the predictor predicted 4, we only use the top 2 predicted oriC to create the matches.' )
    info_lines.append( f'The extra oriCs found by my predictor {did_not} lead to any more matches with oriCs from DoriC.' )
    info_lines.append( f'\tThere were a total of {total_D_oriC} oriCs predicted by DoriC' )
    info_lines.append( f'\tThere were a total of {total_Z_oriC} oriCs predicted by my predictor\n' )

    multi_match_D = [match[0] for match in multi_matching if match[0]]
    multi_match_Z = [match[1] for match in multi_matching if match[1]]
    disp_dist = f'{max_dist*100} % of the sequence length' if max_dist <= 1 else f'{max_dist} bp'

    info_lines.append( f'{total_D_by_Z} out of {total_D_oriC} ({total_D_by_Z/total_D_oriC*100:.2f} %) DoriC oriCs were found by my predictor with max accepted distance of {disp_dist}\n' )

    info_lines.append( f'{len(multi_match_D)} DoriC oriC were matched to mutliple Z_oriC. (> 0 is not okay)' )
    info_lines.append( f'\tThis means there are Z_oriC that are closer than {disp_dist} to each other.') if len(multi_match_D) > 0 else None
    info_lines.append( f'{len(multi_match_Z)} Z_oriC were matched to mutliple DoriC oriC. (> 0 is okay)' )
    info_lines.append( f'\tThis means there are DoriC oriC that are closer than {disp_dist} to each other.\n') if len(multi_match_Z) > 0 else None

    info_lines.append( f'The Z_oriC was off by { int(sum(all_distances["Distance"])/len(all_distances["Distance"])) } bp on average.' )

    if info_file_path is not None:
        with open(info_file_path, 'w') as fh:
            fh.write('\n'.join(info_lines))
    if print_info:
        print('\n'.join(info_lines))

    return pd.DataFrame(all_distances)


if __name__ == '__main__':
    # Input paths
    Z_curve_csv = 'NCBI data prep/oriC_predictions_v2_csvs/All_data.csv'
    DoriC_csv   = 'DoriC data prep/DoriC_oriC_concat_entries.csv'

    # Output paths
    comparator_csv = 'Comparison/v2/in_both_sets_all.csv'
    all_file_path  = 'Comparison/v2/comparison_info_file_all.txt'

    # Make or load csv
    # comparator_df = make_comparator_csv(Z_curve_csv, DoriC_csv, comparator_csv=comparator_csv)
    comparator_df = pd.read_csv(comparator_csv)

    # General info
    avg_seq_len = comparator_df['Sequence_length'].sum() // comparator_df.shape[0]
    disp_dist   = f'{MAX_DIST*100} % of the sequence length' if MAX_DIST <= 1 else f'{MAX_DIST} bp'
    print('Max genome length          :', comparator_df['Sequence_length'].max())
    print('Min genome length          :', comparator_df['Sequence_length'].min())
    print('Mean genome length         :', avg_seq_len)
    print('Deviation in genome length :', int(comparator_df['Sequence_length'].std()))
    print('Genome size is not very consistent. The max distance that two oriC can be from each other will')
    print(f'be set to {disp_dist} of the total genome length rather than a fixed number for all genomes.\n')

    # # ALL DATA
    hist_all_df = compare_dbs(df=comparator_df, info_file_path=all_file_path, print_info=True, max_dist=MAX_DIST)
    # hist_all_df.sort_values(by='Distance', inplace=True)
    # print(hist_all_df.tail(10))
    hist_all_df.to_csv('Comparison/v2/hist_all_df.csv', index=False)
    # # hist_all_df = pd.read_csv('Comparison/hist_all_df.csv')
    # # pf.distance_histogram(hist_all_df, log=True)
    # # pf.distance_histogram(hist_all_df, log=False)

    # # All the same steps, but only the EXPERIMENTAL DATA
    # exp_file_path  = 'Comparison/v2/comparison_info_file_experimental.txt'
    # exp_refseq = [ # Accessions that have been experimentally verified.
    #     'NC_000964', 'NC_002947', 'NC_003272', 'NC_003869', 'NC_003888', 'NC_005090', 'NC_006461', 'NC_007633', 'NC_000913', 'NC_003047',
    #     'NC_007604', 'NC_000962', 'NC_002696', 'NC_002971', 'NC_005363', 'NC_008255', 'NC_009850', 'NC_010546', 'NC_010547', 'NC_011916'
    # ]
    # experiment_df = comparator_df[comparator_df['RefSeq'].isin(exp_refseq)]
    # hist_exp_df = compare_dbs(df=experiment_df, info_file_path=exp_file_path, print_info=True, max_dist=MAX_DIST)
    # # pf.distance_histogram(hist_exp_df, log=True)
    # # pf.distance_histogram(hist_exp_df, log=False)

    # # Only plot those without False_order == False
    # print('False Order\n')
    # order_df = comparator_df[comparator_df['False_order'] == False]
    # hist_order_df = compare_dbs(df=order_df, info_file_path=None, print_info=True, max_dist=MAX_DIST)
    # print()

    # # N-count in genome
    # print('N-Penalty\n')
    # n_list = [x for x in comparator_df['N_penalty']]
    # avg_n = sum(n_list)/len(n_list)
    # hist_n_penalty = compare_dbs(df=comparator_df, info_file_path=None, print_info=True, max_dist=MAX_DIST, n=avg_n)
    # print()

    # print('D-Penalty\n')
    # comparator_df['D_penalty'] = comparator_df['D_penalty'].map(literal_eval)
    # d_list = [y for x in comparator_df['D_penalty'] for y in x]
    # avg_d = sum(d_list)/len(d_list)
    # hist_d_penalty = compare_dbs(df=comparator_df, info_file_path=None, print_info=True, max_dist=MAX_DIST, d=avg_d)
    # print()

    # # Option used to determine oriC
    # print('0-option')
    # zero_hist = compare_dbs(df=comparator_df, info_file_path=None, print_info=True, max_dist=MAX_DIST, alt_oriC='xy_oriC')
    # print()

    # print('1-option')
    # one_df = comparator_df[comparator_df['O_penalty'] == 1]
    # one_hist = compare_dbs(df=one_df, info_file_path=None, print_info=True, max_dist=MAX_DIST)
    # print()

    # print('2-option')
    # two_hist = compare_dbs(df=comparator_df, info_file_path=None, print_info=True, max_dist=MAX_DIST, alt_oriC='xgc_oriC')
    # print()

    # print('3-option')
    # three_hist = compare_dbs(df=comparator_df, info_file_path=None, print_info=True, max_dist=MAX_DIST, alt_oriC='ygc_oriC')
    # print()
