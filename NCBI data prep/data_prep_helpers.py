import csv
import os
import shutil

os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

def merge_csvs(file_folder, merged_csv, fieldnames, length=-1, headers=False):
    '''
    Merge multiple csvs into one csv.
    Arguments:
        file_folder : path to folder with csvs that have to be merged
        merged_csv  : name of the merged csv
        fieldnames  : list of fieldnames for the merged_csv
        length      : amount of rows in each single csv. -1, if the length of each single
                      csv is not known or not the same.
        headers     : if the single csvs have headers or not
    '''
    file_list = os.listdir(file_folder)
    with open(merged_csv, 'w', newline='') as fh_out:
        writer = csv.writer(fh_out)
        writer.writerow(fieldnames)
        for file in file_list:
            with open(file_folder + file, 'r') as fh_in:
                reader = csv.reader(fh_in)
                for i, row in enumerate(reader):
                    if headers and i == 0:
                        pass
                    elif i == length:
                        break
                    else:
                        writer.writerow(row)


def move_fastas(db_loc, on_cluster=True, split=4):
    '''
    Split one folder into 'split' amount of folders with roughly the same amount of files in them.
    Used for easier parallel processing. Instead of one big job with 4000 samples. Run 4 jobs with 1000 samples at the same time.
    '''
    if on_cluster:
        path = db_loc + '/bacteria'
        samples = os.listdir(path)
    else:
        path = db_loc + '/chromosomes_only'
        samples = os.listdir(path)

    samples_per_split = len(samples)//split
    for i in range(split):
        new_folder = db_loc + '/' + str(i)
        os.mkdir(new_folder)
        start = 0 + i*samples_per_split
        stop  = (i+1)*samples_per_split if i != split-1 else None
        for sample in samples[start:stop]:
            shutil.move(path + '/' + sample, new_folder + '/' + sample)


if __name__ == '__main__':

    max_oriCs = 10
    # fieldnames = ['RefSeq', 'Organism', 'Sequence_length', 'N_penalty', 'O_penalty', 'D_penalty', 'False_order']
    # fieldnames.extend([f'oriC_edges_{i}' for i in range(max_oriCs)])
    # fieldnames.extend([f'oriC_middles_{i}' for i in range(max_oriCs)])
    # fieldnames.extend([f'xy_oriC_{i}' for i in range(max_oriCs)])
    # fieldnames.extend([f'xgc_oriC_{i}' for i in range(max_oriCs)])
    # fieldnames.extend([f'ygc_oriC_{i}' for i in range(max_oriCs)])

    fieldnames = ['RefSeq', 'Organism', 'Sequence_length', 'False_order', 'GC_Concentration']
    fieldnames.extend([f'oriC_middles_{i}' for i in range(max_oriCs)])
    fieldnames.extend([f'Occurance_oriC_{i}' for i in range(max_oriCs)])

    merge_csvs('refseq_23/csvs/', 'refseq_23/Predicted_23.csv', fieldnames, length=1)
    # move_fastas('refseq_23', on_cluster=False, split=4)