import re

from Bio import SeqIO, Entrez, SeqRecord
from typing import TextIO, Union

# Self-made module
from peak import Peak

def fetch_file(accession: str, email: str, api_key: Union[str, None], rettype: str) -> TextIO:
    """Downloads the given file_tpye of the given accession in temporary memory"""
    Entrez.email = email
    if api_key is not None:
        Entrez.api_key = api_key
    return Entrez.efetch(db="nuccore", id=accession, rettype=rettype, retmode="text")


def read_FASTA(handle: TextIO) -> SeqRecord.SeqRecord:
    """Read a file and returns the SeqRecord of only the first entry in the file"""
    Seq_records = SeqIO.parse(handle, 'fasta')
    Seq_obj = next(Seq_records)
    return Seq_obj.id, Seq_obj.seq 


def read_gene_info(handle: TextIO , genes_list: list) -> dict:
    """Read FASTA-file acquired with rettype='fasta_cds_na'."""
    obj = SeqIO.parse(handle, 'fasta')
    genes = [gene.lower() for gene in genes_list]
    genes_dict = {}
    for gene in obj:
        features = [x.split('=') for x in re.findall(r"\[(.*?)\]", gene.description) if '=' in x]
        feature_dict = {feature[0] : feature[1] for feature in features}
        try: gene_name = feature_dict.pop('gene')
        except KeyError:
            try: gene_name = feature_dict['protein'].split(' ')[0]
            except KeyError: continue
        if gene_name.lower() in genes:
            genes_dict.update({gene_name : feature_dict})
    return genes_dict


def handle_location(location: str) -> list:
    '''Gene locations come in four flavours, each has to be handled differently.'''
    handled = []
    if 'complement' in location:
        handled.append( location.lstrip('complement(').rstrip(')').split('..') )
    elif 'join' in location:
        locs_list = location.lstrip('join(').rstrip(')').split(',')
        handled.extend( [loc.split('..') for loc in locs_list] )
    elif '<' in location:
        handled.append( location.lstrip('<').split('..') )
    else:
        handled.append( location.split('..') )
    # DEBUG REMOVE LATER
    if len(handled) < 2:
        return None
    return [[int(loc[0]), int(loc[1])] for loc in handled]


def extract_locations(seq_len: int, genes_dict: dict) -> list:
    '''Returns Peaks of the middle position of every gene in the dictionary.'''
    locations = []
    for gene_dict in genes_dict.values():
        locations.extend(handle_location(gene_dict['location']))
    # DEBUG REMOVE LATER
    if None in locations:
        return None
    middles = [Peak(Peak.get_middle(loc[0], loc[1], seq_len), seq_len, 0) for loc in locations]
    return middles