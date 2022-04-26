import ncbi_genome_download as ngd

# https://github.com/kblin/ncbi-genome-download
# See .config file for download parameters. That's the most documentation you'll get for use as a module.
# Don't be alarmed if only MD5SUM5 hashes show up at first. Those get downloaded first, the fna.gz files after.

ngd.download(section='refseq', file_formats='fasta', refseq_categories='reference', assembly_levels='complete', groups='bacteria')
# refseq_categories='reference' downloads 15 sequences (24-04-2022)
# refseq_categories='representative' downloads ~15k sequences (24-04-2022) (I didn't finish the download)