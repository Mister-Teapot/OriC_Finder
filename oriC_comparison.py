import os
from ast import literal_eval
import pandas as pd

# Set cwd to location of this script
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

# Pandas printing options
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 15)

def compare_dbs(Z_curve_csv, DoriC_csv, comparator_csv=None):

    # Read CSVs
    Z_curve_df = pd.read_csv( Z_curve_csv )
    DoriC_df   = pd.read_csv( DoriC_csv )

    # Completely useless! Dummy...
    # Setting columns to proper types: Z-curve
    # Z_curve_df = Z_curve_df.astype({'Sequence_length':'float', 'False_order':'bool'})
    # Z_curve_df['Penalties'] = Z_curve_df['Penalties'].map(literal_eval)

    # oriC_cols_Z = [i for i in Z_curve_df.columns if 'oriC' in i]
    # most_oriC_Z = int(sorted(oriC_cols_Z, key=lambda x:x[-1], reverse=True)[0][-1]) + 1
    # for i in range(most_oriC_Z):
    #     Z_curve_df = Z_curve_df.astype({f'oriC_middles_{i}':'float'})
    #     Z_curve_df[f'oriC_edges_{i}'].fillna('()', inplace=True)
    #     Z_curve_df[f'oriC_edges_{i}'] = Z_curve_df[f'oriC_edges_{i}'].map(literal_eval)

    # # Setting columns to proper types: DoriC
    # oriC_cols_D = [i for i in DoriC_df.columns if 'oriC' in i]
    # most_oriC_D = int(sorted(oriC_cols_D, key=lambda x:x[-1], reverse=True)[0][-1]) + 1
    # for i in range(most_oriC_D):
    #     DoriC_df[f'DoriC_oriC_{i}'].fillna('()', inplace=True)
    #     DoriC_df[f'DoriC_oriC_{i}'] = DoriC_df[f'DoriC_oriC_{i}'].map(literal_eval)

    # Merge into single dataframe
    col_names_Z = [x for x in Z_curve_df.columns if x not in ['RefSeq', 'Organism']]
    in_both_df = {x: [] for x in DoriC_df.columns.to_list() + col_names_Z}
    for i, sample in Z_curve_df.iterrows():
        if sample['RefSeq'] in DoriC_df['RefSeq'].to_list():
            for j in col_names_Z:
                in_both_df[j].append( sample[j] )
            for j in DoriC_df.columns.to_list():
                in_both_df[j].append( DoriC_df.loc[DoriC_df['RefSeq'] == sample['RefSeq'], j].values[0] ) # Series.values -> ['element'], very annoying

    del Z_curve_df
    del DoriC_df

    comparator_df = pd.DataFrame(in_both_df).head()
    if comparator_csv is not None:
        comparator_df.to_csv(comparator_csv, index=False)

    # Compare predictions (Assuming DoriC as ground truth for now)
    num_of_predictions = [0, 0] # [more oriC in DoriC, more oriC by mine]
    for i, sample_original in comparator_df.iterrows():
        sample = sample_original.copy()
        sample.dropna(inplace=True)

        num_oriC_D = int( sorted([i for i in sample.axes[0] if 'DoriC_oriC' in i], key=lambda x:x[-1], reverse=True)[0][-1] ) + 1
        num_oriC_Z = int( sorted([i for i in sample.axes[0] if 'oriC_edges' in i], key=lambda x:x[-1], reverse=True)[0][-1] ) + 1

        if num_oriC_D > num_oriC_Z:
            num_of_predictions[0] += 1
            ...
        elif num_oriC_D < num_oriC_Z:
            num_of_predictions[1] += 1
            ...
        else: # Equal amount of oriC predictions
            ...

if __name__ == '__main__':
    Z_curve_csv = 'NCBI data prep/refseq_15/NCBI_oriC_15_improved.csv'
    DoriC_csv   = 'DoriC data prep/DoriC_oriC_concat_entries.csv'
    compare_dbs(Z_curve_csv, DoriC_csv, comparator_csv='Both_predictions.csv')