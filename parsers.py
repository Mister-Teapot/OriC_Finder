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


def read_gene_info(handle: TextIO, genes_list: list) -> dict:
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
    '''Gene locations come in a billion flavours due to annotation options. This function handles >95 % of cases.'''
    handled = []
    loc = ''.join([x for x in location if x not in ['<', '>']]) # remove fuzziness from location
    raw_locs = loc.split(',')
    locs = []
    print(raw_locs)
    done, i = False, 0
    while not done:
        # Not perfect yet, assumes join-group only has two entries
        # if element after the start of a group has balanced brackets, then sth is up.
        if i == len(raw_locs) - 1:
            done = True
        if raw_locs[i].count('(') != raw_locs[i].count(')'):
            locs.append(raw_locs[i] + ',' + raw_locs[i+1])
            i += 2
        else:
            locs.append(raw_locs[i])
            i += 1
    print(locs)

    '''
    raw_entries = re.split(',\s?', location) # dont even trust them to comma-separate without spaces
    for entry in raw_entries:
        raw_loc = re.search(r'\d+(\.\.[<>]?\d*)?', entry).group() # SHOULD only be one group
        loc_coords = re.split(r'\.\.[<>]?', raw_loc)
        if len(loc_coords) == 1:
            loc_coords.append(loc_coords[0])
        elif len(loc_coords) == 0:
            print('Something went horribly wrong. Gene has been very weirdly annotated or gene did not get a location\n', location)
        # ONE MORE safety check for any non-numerical characters in the location
        for i in range(len(loc_coords)):
            loc_coords[i] = re.search(r'\d+', loc_coords[i]).group() # should DEFINITELY be only one group
        handled.append(loc_coords)
    try:
        to_return = [[int(loc[0]), int(loc[1])] for loc in handled]
    except: # If this becomes something other than `return None`, then implement separate KeyboardInterrupt
        return None
    return to_return
    '''

def extract_locations(seq_len: int, genes_dict: dict) -> list:
    '''Returns Peaks of the middle position of every gene in the dictionary.'''
    locations = []
    for gene_dict in genes_dict.values():
        clean_locs = handle_location(gene_dict['location'])
        if clean_locs is None:
            return None
        locations.extend(clean_locs)
    middles = [Peak(Peak.get_middle(loc[0], loc[1], seq_len), seq_len, 0) for loc in locations]
    return middles