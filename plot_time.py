import sys
import pandas as pd
import multiprocessing as mp
import matplotlib.pyplot as plt

sys.path.append('../OriC_Finder/')

from oriC_Finder import find_oriCs

def make_comparator_csv():
    df = pd.read_csv('NCBI data prep/oriC_predictions_v5_csvs/Predicted_DoriC_accessions.csv')
    df.sort_values(by='Sequence_length', inplace=True)
    spacing = [i for i in range(0, df.shape[0], df.shape[0]//100)] + [df.shape[0]-1]

    accessions = df.iloc[spacing]['RefSeq'].tolist()

    plot_dict = {
        'seq_len': [],
        'num_of_genes': [],
        'calc_disp_time': [],
        'read_genes_time': [],
        'total_time': []
    }

    for i, acc in enumerate(accessions):
        properties = find_oriCs(accession=acc, email='zoyavanmeel@gmail.com')

        plot_dict['seq_len'].append(properties['seq_size'])
        plot_dict['num_of_genes'].append(properties['num_of_genes'])
        plot_dict['calc_disp_time'].append(properties['time_of_prediction'][0])
        plot_dict['read_genes_time'].append(properties['time_of_prediction'][1])
        plot_dict['total_time'].append(properties['time_of_prediction'][2])
        print(f'{i+1}/102 done: {acc}')

    pd.DataFrame(plot_dict).to_csv('Time_performance.csv')

df = pd.read_csv('Time_performance.csv')
plt.scatter(df['seq_len'], df['total_time'])
plt.show()
plt.scatter(df['num_of_genes'], df['total_time'])
plt.show()
