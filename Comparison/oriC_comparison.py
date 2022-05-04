import os, sys
from ast import literal_eval
from itertools import product
import pandas as pd

# Pandas printing options
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 15)

# Self-made module
#   Local path
sys.path.append('../OriC_Finder')
#   Cluster path
# sys.path.append('/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/OriC_Finder/')
from oriC_Finder import get_peak_windows, merge_peaks


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

def compare_dbs(df=None, csv=None, info_file_path=None, print_info=False):
    '''
    Compares results. Results can be loaded from a pandas.DataFrame OR csv.
    info       : if str, writes a text_file to the given location with some info on the analysis.
    print_info : if True, prints the same info that was written to 'info' with print-statements for quick reading.
    '''

    if df is None and csv is None:
        raise KeyError('Must load a pandas.DataFrame or CSV file.')
    if df is not None and csv is not None:
        raise KeyError('Load either a pandas.DataFrame or CSV, not both.')

    if csv is not None:
        df = pd.read_csv(csv)

    # Z_curve_df['Penalties'] = Z_curve_df['Penalties'].map(literal_eval) # For converting tuple-strings to tuple objects

    # Compare predictions (Assuming DoriC as ground truth for now)
    num_of_predictions = [0, 0] # [more oriC in DoriC, more oriC by mine]
    for i, sample_original in df.iterrows():
        sample = sample_original.copy()
        sample.dropna(inplace=True)

        D_oriC_cols = sorted( [i for i in sample.axes[0] if 'DoriC_oriC' in i], key=lambda x:x[-1] )
        Z_oriC_cols = sorted( [i for i in sample.axes[0] if 'oriC_edges' in i], key=lambda x:x[-1] )

        if len(D_oriC_cols) > len(Z_oriC_cols):
            num_of_predictions[0] += 1
        elif len(D_oriC_cols) < len(Z_oriC_cols):
            num_of_predictions[1] += 1

        for oriC_D_col, oriC_Z_col in product(D_oriC_cols, Z_oriC_cols):
            print(literal_eval(sample[oriC_D_col]))

    info_lines = []
    info_lines.append( f'There were {df.shape[0]} organisms.' )
    info_lines.append( f'\t{df.shape[0] - sum(num_of_predictions)} times both databases predicted the same amount of oriCs.' )
    info_lines.append( f'\t{num_of_predictions[0]} times DoriC predicted more oriCs.' )
    info_lines.append( f'\t{num_of_predictions[1]} times mine predicted more oriCs.' )

    if info_file_path is not None:
        with open(info_file_path, 'w') as fh:
            fh.write('\n'.join(info_lines))
    if print_info:
        print('\n'.join(info_lines))

if __name__ == '__main__':
    Z_curve_csv    = 'NCBI data prep/299+15_merged.csv'
    DoriC_csv      = 'DoriC data prep/DoriC_oriC_concat_entries.csv'

    comparator_csv = 'Comparison/Both_predictions_out_of_299+15.csv'
    info_file_path = 'Comparison/test_info_file.txt'

    comparator_df = make_comparator_csv(Z_curve_csv, DoriC_csv, comparator_csv=comparator_csv)
    compare_dbs(df=comparator_df, info_file_path=info_file_path, print_info=True) # csv=comparator_csv