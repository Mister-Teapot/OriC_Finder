import ncbi_genome_download as ngd

'''
https://github.com/kblin/ncbi-genome-download
See .config file for download parameters. That's the most documentation you'll get for use as a module.
Don't be alarmed if only MD5SUM5 hashes show up at first. Those get downloaded first, the fna.gz files after.

parallel=2 for parallel downloading (set up multithread in .sbatch)

refseq_categories = [
    'reference',
    'representative',
    'na'
]
'na': most are not in either other category. Have been entirely sequenced, but are not used as a reference and do not represent 

dry_run=True to get the amount of genomes your settings would download.
'''

# ngd.download(section='refseq', file_formats='fasta', groups='bacteria', assembly_levels='complete', dry_run=True) # 26516
# ngd.download(section='refseq', file_formats='fasta', groups='bacteria', assembly_levels='complete', refseq_categories='na', dry_run=True) # 22930
# ngd.download(section='refseq', file_formats='fasta', groups='bacteria', assembly_levels='complete', refseq_categories='representative', dry_run=True) # 3571
# ngd.download(section='refseq', file_formats='fasta', groups='bacteria', assembly_levels='chromosome', dry_run=False) # 4222

# scaffolds and contigs are not eligible for this projecct.
# ngd.download(section='refseq', file_formats='fasta', groups='bacteria', assembly_levels='scaffold', dry_run=True) # 82496
# ngd.download(section='refseq', file_formats='fasta', groups='bacteria', assembly_levels='contig', dry_run=True) # 132482