'''
Dataset structure via website:
genome_assemblies_genome_fasta/                       # 1. You actually get a .tar-file. Unzip that first
├── ncbi-genomes-YYYY-MM-DD/                          # 2. folder with datestamp
│   ├── GCF_000006925.2.fna.gz                        # 3. zip file that contains the .fasta-file with the genome sequence
│   │   └── GCF_000006925.2_ASM692v2_genomic.fna      # 4. the actual fasta file we want...
│   ├── ...
│   ├── md5checksums.txt                              # -. Collection of MD5SUM5 hashes
│   └── README.txt                                    # -. Generic readme, explains download options
└── report.txt

Dataset structure via FTP:
refseq/                                               # 1. the folder you get
├── bacteria/                                         # 2. layer of all groups you chose to download
│   ├── GCF_000006925.2/                              # 3. folders with names of all the species of said group you chose to download
│   │   ├── GCF_000006925.2_ASM692v2_genomic.fna.gz   # 4. zip file that contains the .fasta-file with the genome sequence
│   │   │   └── GCF_000006925.2_ASM692v2_genomic.fna  # 5. the actual fasta file we want...
│   │   └── MD5SUM5                                   # -. MD5SUM5 hash file
|   └── .../


After running this script the database will be structured as follows, unless a different out_loc was specified:
refseq/
├── bacteria/
│   ├── GCF_000006925.2_ASM692v2_genomic.fna
│   └── ...

'''

import gzip
import os
import shutil

os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

def read_database(db_loc, method, inplace=True, split_genome=False, out_loc=None, filter_plasmids=False):
    """
    Extracts FASTA files from the raw NCBI download. Only works if you only downloaded FASTA files.
    Arguments:
        db_loc          : path to location of database
        method          : ['website', 'FTP']
                          The structure of the download you get directly from the website is different
                          than that of the one you get with the FTP (eg. ngd module, datasets commandline tool)
                              'website' : small dataset downloaded directly from the NCBI website
                              'FTP'     : downloaded via the NCBI FTP
        inplace         : if True, deletes the original download to save space.
                          NOTE: PERMANENTLY deletes the original folder with the gz and hash files. Set to False if memory
                          space is not an issue. I added it for my personal benefit otherwise I would need to temporarily
                          store the database twice (not always possible) and delete the compressed folders manually
        split_genome    : if True, gives each chromosome and plasmid a separate FASTA file.
                          NOTE: Files that need splitting will get new names AND happens in the place of the out_loc !!
                          New name = old name + _i, with i = position in original FASTA.
        filter_plasmids : if True, deletes all FASTA files with plasmid sequences.
                          NOTE: split_genome must be True for filter to work.
        out_loc         : Extracted dataset output path. Necessary if inplace=False.
    """

    method_options = [
        'website', 
        'FTP'
    ]


    # Some error handling
    if method not in method_options:
        raise ValueError('Method not in method types.')
    if inplace == False and out_loc is None:
        raise ValueError('No output path given')
    if filter_plasmids and not split_genome:
        raise ValueError('Cannot remove plasmids without splitting genomes into seperate files.')


    if method == 'website':
        # Unzipping this is very similar to FTP
        db = os.path.join(db_loc, 'genome_assemblies_genome_fasta')
        if out_loc is None:
            os.mkdir( os.path.join( db, 'bacteria') )
            out_loc = os.path.join( db, 'bacteria')

        datestamp = os.listdir(db)[1]           # ['bacteria', 'ncbi-genomes-YYYY-MM-DD', 'report.txt']
        path = os.path.join(db, datestamp)
        sample_names = os.listdir(path).copy()  # ['accession_num.fna.gz', ...]
        sample_names.remove('md5checksums.txt')
        sample_names.remove('README.txt')
        for sample in sample_names:
            fna_file = sample[:-3]              # 'accession_num.fna'
            with gzip.open( os.path.join(path, sample), 'rb' ) as f_in, open( os.path.join(out_loc, fna_file), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        # !!! PERMANENTLY deletes the original folder with the gz and hash files. Remove if memory space is not an issue.
        # I added it for my personal benefit otherwise I would need to temporarily store the database twice and delete the compressed folders manually
            if inplace:
                os.remove(os.path.join(path, sample))
        if inplace:
            shutil.rmtree(path)
        os.rename(db, os.path.join(db_loc, 'refseq'))


    if method == 'FTP':
        db_loc = os.path.join(db_loc, 'refseq', 'bacteria')
        if out_loc is None:
            out_loc = db_loc
        sample_names = os.listdir(db_loc).copy()
        for sample in sample_names:
            path     = os.path.join(db_loc, sample)
            gz_file  = os.listdir(path)[0]              # 'accession_num.fna.gz'
            fna_file = gz_file[:-3]                     # 'accession_num.fna'

            with gzip.open( os.path.join(path, gz_file), 'rb' ) as f_in, open( os.path.join(out_loc, fna_file), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
            # !!! PERMANENTLY deletes the original folder with the gz and hash files. Remove if memory space is not an issue.
            # I added it for my personal benefit otherwise I would need to temporarily store the database twice and delete the compressed folders manually
            if inplace:
                shutil.rmtree(path)


    if split_genome:
        # Read through the FASTA twice. Once to assert how many new files are needed, once to place the right info in each file.
        # There is probably a much better way to do this, but I don't know how to safely do it. While-loops for file handling sound scary.
        in_loc = out_loc
        sample_list = os.listdir(in_loc).copy()
        for sample in sample_list:
            files_needed = 0
            with open( os.path.join(in_loc, sample), 'r' ) as fh:
                for line in fh:
                    if line[0] == '>':
                        files_needed += 1

            start_line = 0
            for i in range(files_needed):
                with open( os.path.join(in_loc, sample), 'r' ) as f_in, open( os.path.join(out_loc, sample[:-4]+'_'+str(i)+sample[-4:]), 'w' ) as f_out:
                    name_found = False
                    for j, line in enumerate(f_in):
                        if j >= start_line:
                            if name_found and line[0] == '>':
                                start_line = j
                                break
                            elif line[0] == '>':
                                f_out.write(line)
                                name_found = True
                            else:
                                f_out.write(line)
            os.remove(os.path.join(in_loc, sample))


    if filter_plasmids:
        in_loc = out_loc
        sample_list = os.listdir(in_loc).copy()
        for sample in sample_list:
            delete = False
            with open( os.path.join(in_loc, sample), 'r' ) as fh:
                for line in fh:
                    if 'plasmid' in line:
                        delete = True
            if delete:
                os.remove(os.path.join(in_loc, sample))


if __name__ == "__main__":
    read_database(r"C:\0. School\Bachelor Thesis\Zoya_Code+Data\oriFinder\oriFinder Comparison\NCBI_data_prep", method='FTP', split_genome=True, inplace=True, filter_plasmids=True)