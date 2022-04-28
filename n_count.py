import os
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

def read_FASTA(filename):
    """Read in a file in FASTA format and returns the name and sequence"""
    n_count = 0
    name = None
    with open(filename, 'r') as fh:
        for line in fh:
            line = line.rstrip()
            if line[0] == '>' and name is None:
                name = line[1:]
            elif line[0] == '>' and name is not None:
                break # In case the FASTA contains multiple sequences, only read the first one.
            else:
                line = line.upper()
                n_count += line.count('N')
    return name, n_count

db_loc = './NCBI data prep//refseq_287/chromosomes_only'
count = []
for fasta in os.listdir(db_loc):
    _, single_count = read_FASTA(os.path.join(db_loc, fasta))
    count.append(single_count)
print('Average N-count in 287 representative refseq genomes: ', sum(count)/len(count))