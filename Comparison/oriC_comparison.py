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

MAX_DIST = 1/50


def make_comparator_csv(Z_curve_csv, DoriC_csv, comparator_csv):
    '''Merges Z_curve_csv and DoriC_csv into one csv and/or dataframe for comparison with compare_dbs'''

    # Read CSVs
    Z_curve_df = pd.read_csv( Z_curve_csv )
    DoriC_df   = pd.read_csv( DoriC_csv )

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

def compare_dbs(df=None, csv=None, info_file_path=None, print_info=False, max_dist=30000):
    '''
    Compares results. Results can be loaded from a pandas.DataFrame OR csv.
    info       : If str, writes a text_file to the given location with some info on the analysis.
    print_info : If True, prints the same info that was written to 'info' with print-statements for quick reading.
    max_dist   : Maximum distance two oriCs can have from each other if they want to be matched. If <= 1, will be
                 seen as a fraction of the total sequence length. If > 1, it will be seen as # of basepairs.
    '''

    if df is None and csv is None:
        raise KeyError('Must load a pandas.DataFrame or CSV file.')
    if df is not None and csv is not None:
        raise KeyError('Load either a pandas.DataFrame or CSV, not both.')

    if csv is not None:
        df = pd.read_csv(csv)
    
    if max_dist > 30000:
        warnings.warn("Warning: Raising the max_dist above the minimum distance of oriC found with the predictor. Potential overlap in DoriC to Z_oriC matching")

    # Z_curve_df['Penalties'] = Z_curve_df['Penalties'].map(literal_eval) # For converting tuple-strings to tuple objects

    # Compare predictions (Assuming DoriC as ground truth for now)
    num_of_predictions = [0, 0] # [more oriC in DoriC, more oriC by mine]
    multi_matching = []
    list_all_distances = []
    total_D_oriC = 0
    total_Z_oriC = 0
    total_D_by_Z = 0
    for i, sample_original in df.iterrows():
        sample = sample_original.copy()
        sample.dropna(inplace=True)

        seq_len = sample['Sequence_length']
        D_oriC_cols = sorted( [i for i in sample.axes[0] if 'DoriC_oriC' in i], key=lambda x:x[-1] )
        Z_oriC_cols = sorted( [i for i in sample.axes[0] if 'oriC_middle' in i], key=lambda x:x[-1] )
    
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

        all_distances = [] # D_oriC index and its closest Z_oriC
        for i, D_oriC in enumerate(D_oriC_middles):
            distances = []
            for j, Z_oriC in enumerate(Z_oriC_middles[:len(D_oriC_middles)]):
                # Modified: always take D_oriC as 0-point
                dist_1 = D_oriC - Z_oriC
                dist_2 = Z_oriC + seq_len-1 - D_oriC
                distances.append( (j, min(dist_1, dist_2, key=lambda x:abs(x))) )
                # Original: always positive distance
                # dist_1 = max(D_oriC, Z_oriC) - min(D_oriC, Z_oriC)
                # dist_2 = min(D_oriC, Z_oriC) + seq_len-1 - max(D_oriC, Z_oriC)
                # distances.append( (j, min(dist_1, dist_2)) )
            all_distances.append( (i, min(distances, key=lambda x:x[1])) )
        list_all_distances.append([(x[0], x[1][0], x[1][1]) for x in all_distances])


        ### Accuracy per sample: how many samples were properly matched:

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
    info_lines.append( f'\tThis means there are DoriC oriC that are closer than {disp_dist} to each other.') if len(multi_match_Z) > 0 else None


    if info_file_path is not None:
        with open(info_file_path, 'w') as fh:
            fh.write('\n'.join(info_lines))
    if print_info:
        print('\n'.join(info_lines))


if __name__ == '__main__':
    Z_curve_csv    = 'NCBI data prep/3k+299+15_merged.csv'
    DoriC_csv      = 'DoriC data prep/DoriC_oriC_concat_entries.csv'

    comparator_csv = 'Comparison/Both_predictions_out_of_3k+299+15.csv'
    info_file_path = 'Comparison/comparison_info_file_3k+299+15.txt'
    exp_file_path  = 'Comparison/comparison_info_file_experimental.txt'

    # comparator_df = make_comparator_csv(Z_curve_csv, DoriC_csv, comparator_csv=comparator_csv)
    comparator_df = pd.read_csv('Comparison/Both_predictions_out_of_3k+299+15.csv')
    avg_seq_len = comparator_df['Sequence_length'].sum() // comparator_df.shape[0]

    print('Max genome length          :', comparator_df['Sequence_length'].max())
    print('Min genome length          :', comparator_df['Sequence_length'].min())
    print('Mean genome length         :', avg_seq_len)
    print('Deviation in genome length :', int(comparator_df['Sequence_length'].std()))
    print('Genome size is not very consistent. The max distance that two oriC can')
    disp_dist = f'{MAX_DIST*100} % of the sequence length' if MAX_DIST <= 1 else f'{MAX_DIST} bp'
    print(f'be from each other will be set to {disp_dist} of the total genome length rather than a fixed number for all genomes.\n')

    # compare_dbs(df=comparator_df, info_file_path=info_file_path, print_info=True, max_dist=MAX_DIST)

    # Accessions that have been experimentally verified.
    exp_refseq = [
        'NC_000964', 'NC_002947',
        'NC_003272', 'NC_003869',
        'NC_003888', 'NC_005090',
        'NC_006461', 'NC_007633',
        'NC_000913', 'NC_003047',
        'NC_007604', 'NC_000962',
        'NC_002696', 'NC_002971',
        'NC_005363', 'NC_008255',
        'NC_009850', 'NC_010546',
        'NC_010547', 'NC_011916'
    ]

    print(comparator_df[comparator_df['RefSeq'].isin(exp_refseq)]['RefSeq'])
    exp_df = comparator_df[comparator_df['RefSeq'].isin(exp_refseq)]

    # compare_dbs(df=exp_df, info_file_path=exp_file_path, print_info=True, max_dist=MAX_DIST)
