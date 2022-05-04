import os
from ast import literal_eval
from itertools import product
import pandas as pd

# Set cwd to location of this script
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

# Pandas printing options
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 15)

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

    comparator_df = pd.DataFrame(in_both_df).head()
    comparator_df.to_csv(comparator_csv, index=False)

    return comparator_df

def compare_dbs(df=None, csv=None):
    '''Compares results. Results can be loaded from a pandas.DataFrame or csv.'''

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
        for oriC_D, oriC_Z in product(D_oriC_cols, Z_oriC_cols):
            ...

if __name__ == '__main__':
    Z_curve_csv    = 'NCBI data prep/refseq_15/NCBI_oriC_15_improved.csv'
    DoriC_csv      = 'DoriC data prep/DoriC_oriC_concat_entries.csv'
    comparator_csv = 'Both_predictions.csv'

    comparator_df = make_comparator_csv(Z_curve_csv, DoriC_csv, comparator_csv=comparator_csv)
    compare_dbs(df=comparator_df) # csv=comparator_csv