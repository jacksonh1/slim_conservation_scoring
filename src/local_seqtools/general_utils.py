import copy
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from Bio import Align, AlignIO, Seq, SeqIO

sys.path.append("/Users/jackson/tools/iupred2a")
sys.path.append("/home/jackson/tools/iupred2a")
sys.path.append("/home/jch/tools/iupred2a")
import iupred2a_lib as iup

# from sklearn import metrics


def parse_filename(filename, regex = r'(?P<species_id>\d+_\d+)_(?P<odbid>[\w\d]+)_(?P<level>\w+)_(?P<og_id>\d+at\d+)'):
    p = re.compile(regex)
    m = p.match(filename)
    odb_gene_id = m.group('species_id') + ':' + m.group('odbid')
    level = m.group('level')
    og_id = m.group('og_id')
    return odb_gene_id, level, og_id

# ==============================================================================
# // sequence tools
# ==============================================================================
def str2seqrecord(seq_str, id=None, description=None):
    """convert a string to a SeqRecord object

    Parameters
    ----------
    seq_str : str
        sequence string
    id : str, optional
        id of the SeqRecord object, by default None
    description : str, optional
        description of the SeqRecord object, by default None

    Returns
    -------
    SeqRecord
        SeqRecord object with the sequence string
    """    
    return SeqIO.SeqRecord(Seq.Seq(seq_str), id=id, description=description)


def slice_seq(seq:str, st:int, end:int, flank:int=0) -> str:
    """
    slice a sequence at the specified positions
    flank the slice with the specified number of residues to the left and right
    if the flank would go out of bounds, slice to the end/beginning of the sequence and pad with dashes
    
    Parameters
    ----------
    seq : str
        sequence to be sliced
    st : int
        start position of the slice
    end : int
        end position of the slice
    flank : int, optional
        number of residues to flank the slice with, by default 0
    
    Returns
    -------
    str
        sliced sequence
    
    Examples
    --------
    >>> slice_seq('ABCDEF', 2, 4, flank=1)
    >>> s = '012FPpPp890'
    >>> st=10
    >>> end=7
    >>> flank=0
    >>> print(s[st:end+1])
    >>> for c, i in enumerate(slice_seq(s, st, end+1, flank)):
            print(c, i)
    >>> print(slice_seq(s, st, end+1, flank))

    """
    assert st <= end, 'start must be less than or equal to end'
    assert flank >= 0, 'flank must be positive'

    if st-flank < 0:
        left = '-'*abs(flank-st) + seq[:max(0, st)]
    else:
        left = seq[st-flank:st]
    if end+flank > len(seq):
        right = seq[max(0, st):len(seq)] + '-'*(end+flank-len(seq))
    else:
        right = seq[max(0, st):end+flank]
    # print(left)
    # print(right)
    return left+right
    

def gen_kmers(seq, k):
    '''
    generates list of length k "kmers" comprising `seq`

    Parameters
    ----------
    seq : str
        input sequence to split into k-mers. Should be entry in fasta file
    k : int
        length of fragments (k-mers)
    
    Returns
    -------
    kmers : list
        list of strings - k-mers generated from seq
    '''
    k2 = k - 1
    kmers=[]
    for i in range(len(seq) - k2):
        kmers.append(seq[i:i + k])
    return kmers


def import_aacons_scores(aacons_file):
    scores_dict = {}
    with open(aacons_file) as f:
        for line in f:
            linesplit = line.split()
            if len(linesplit) < 1:
                continue
            scores_dict[linesplit[0].strip('#')] = [float(i) for i in linesplit[1:]]
    return scores_dict


def join_overlapping_paths(p1, p2):
    p1 = str(p1)
    p2 = str(p2)
    return p1 + "/" + "/".join([i for i in p2.split("/") if i not in p1.split("/")])


def import_fasta(fasta_path, output_format='list') -> list[SeqIO.SeqRecord]|dict[str, SeqIO.SeqRecord]:
    """import fasta file into a list or dictionary of SeqRecord objects

    Parameters
    ----------
    fasta_path : str
        file path to fasta file
    output_format : str, optional
        output_format of output. Either 'list' or 'dict'. by default 'list'

    Returns
    -------
    list or dictionary
        list or dictionary of SeqRecord objects for each sequence in the fasta file
    """
    allowed_formats = ['list', 'dict']
    with open(fasta_path) as handle:
        if output_format == 'list':
            seqs = list(SeqIO.parse(handle, 'fasta'))
        elif output_format == 'dict':
            seqs = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        else:
            raise ValueError(f"Invalid output format - {output_format}. Expected one of: {allowed_formats}")
    return seqs


class FastaImporter:
    """import fasta file and return seqrecord objects in various formats

    Parameters
    ----------
    fasta_path : str
        file path to fasta file
    """
    def __init__(self, fasta_path: str|Path):
        self.fasta_path = fasta_path

    def import_as_list(self) -> list[SeqIO.SeqRecord]:
        """return list of SeqRecord objects for each sequence in the fasta file

        Returns
        -------
        List[SeqIO.SeqRecord]
            list of SeqRecord objects
        """        
        with open(self.fasta_path) as handle:
            return list(SeqIO.parse(handle, 'fasta'))

    def import_as_dict(self) -> dict[str, SeqIO.SeqRecord]:
        """return dictionary of SeqRecord objects for each sequence in the fasta file

        Returns
        -------
        dict[str, SeqIO.SeqRecord]
            dictionary of SeqRecord objects, keys are the sequence ids and values are the SeqRecord objects
        """        
        with open(self.fasta_path) as handle:
            return SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))


def strip_dashes_from_sequences(sequences):
    """remove `-` characters from sequences in a list of SeqRecord objects

    Parameters
    ----------
    sequences : list
        list of SeqRecord objects

    Returns
    -------
    list
        list of SeqRecord objects with `-` characters removed from sequences
    """    
    sequences = copy.deepcopy(sequences)
    for s in sequences:
        s.seq = Seq.Seq(str(s.seq).replace('-', ''))
    return sequences


def write_fasta(sequences, output_filename):
    """save list of SeqRecord objects to a fasta file

    Parameters
    ----------
    sequences : list
        list of SeqRecord objects
    output_filename : str
        file path to new output fasta file
    """    
    with open(output_filename, 'w') as handle:
        SeqIO.write(sequences, handle, 'fasta')


def longest_common_substring(s1, s2):
    '''find the largest substring shared by two strings'''
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]


# ==============================================================================
# // alignment tools
# ==============================================================================

def reindex_alignment(record):
    """convert indexes of an ungapped sequence to indexes of it's sequence with gaps

    Parameters
    ----------
    record : SeqRecord object
        biopython sequence record object with gaps present

    Returns
    -------
    str
        ungapped sequence
    list
        list of indexes of nongap characters in the gapped sequence. index the return list with 
        non-gap indexes to get the index of the gapped sequence. If the gapped sequence is `A--A--A`, the return list will be `[0, 3, 6]`
    """    
    s = str(record.seq)
    unal_seq = ''
    index_map = []
    for al_pos,i in enumerate(s):
        # print(al_pos, i)
        if i != '-':
            unal_seq = unal_seq + i
            index_map.append(al_pos)
    return unal_seq, index_map


def reindex_alignment_str(seq_str):
    """convert indexes of an ungapped sequence to indexes of it's sequence with gaps

    Parameters
    ----------
    seq_str : str
        string of a sequence with gaps present

    Returns
    -------
    str
        ungapped sequence
    list
        list of indexes of nongap characters in the gapped sequence. index the return list with 
        non-gap indexes to get the index of the gapped sequence. If the gapped sequence is `A--A--A`, the return list will be `[0, 3, 6]`
    """    
    unal_seq = ''
    index_map = []
    for al_pos,i in enumerate(seq_str):
        # print(al_pos, i)
        if i != '-':
            unal_seq = unal_seq + i
            index_map.append(al_pos)
    return unal_seq, index_map


def find_subseq_in_aligned_seq_str(unaligned_subseq_str, aligned_seq_str):
    """
    Find the position of unaligned_subsequence in aligned_seq_str
    return the positions of the subsequence in the aligned sequence indeces
    returns:
        al_start: start position of the subsequence in the aligned sequence
        al_end: end position of the subsequence in the aligned sequence
        aligned_seq_str[al_start:al_end]: the subsequence in the aligned sequence
    """
    unaligned_seq_str, alignment_index_key = reindex_alignment_str(aligned_seq_str)
    if unaligned_subseq_str not in unaligned_seq_str:
        raise ValueError(f"{unaligned_subseq_str} not in `aligned_seq_str`")
    unal_start = unaligned_seq_str.find(unaligned_subseq_str)
    unal_end = unal_start + len(unaligned_subseq_str) - 1
    al_start = alignment_index_key[unal_start]
    al_end = alignment_index_key[unal_end]
    return al_start, al_end, aligned_seq_str[al_start:al_end+1]



# ==============================================================================
# // iupred
# ==============================================================================

def iupred_scores_from_seqrecord(seqrecord):
    return iup.iupred(str(seqrecord.seq))[0]


def iupred_scores_from_seqstring(seqstring):
    return iup.iupred(seqstring)[0]


def get_idr_indexes_untouched_indexing(score_list, cutoff):
    idr_indexes = []
    in_idr = False
    st, end = 0,0
    for c,i in enumerate(score_list):
        # print(c,i)
        if i>cutoff:
            if in_idr:
                end = c
            else:
                st = c
                end = c
                in_idr = True
        else:
            if in_idr:
                if st != end:
                    idr_indexes.append([st,end])
                in_idr = False
            else:
                # print(f'{c}-pass')
                pass
    if in_idr:
        idr_indexes.append([st,end])
    return idr_indexes


# ==============================================================================
# // parse accession ids
# ==============================================================================

def split_uniprot(prot_id):
    j = re.compile(r"^[st][pr]\|(.+)\|(.+)")
    prot = j.findall(prot_id)[0]
    accession= prot[0]
    name=prot[1]
    return name, accession


def get_description(desc):
    r = re.compile(r"^sp\|.+\|\w+ (.+) OS=.+")
    return r.findall(desc)[0]


def convert_id_to_name(uniprot_id):
    import requests
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
    response = requests.get(url)

    if response.status_code == 200:
        content = response.text
        start = content.find("<name>") + len("<name>")
        end = content.find("</name>")
        uniprot_name = content[start:end].strip()
        return uniprot_name
    else:
        return "not_found"


# ==============================================================================
# // regex utilities
# ==============================================================================

def regex2overlapping(regex):
    new_regex = r"(?=(" + regex + "))"
    return new_regex

def get_regex_matches(regex_pattern: str, seq_str: str):
    """searches for all matches of a regex pattern in a sequence string
    returns a generator object that yields the match sequence, start index, and end index

    Parameters
    ----------
    regex_pattern : str
        regular expression pattern
    seq_str : str
        string to search for matches

    Yields
    ------
    tuple
        (match sequence, start index, end index)
    """    
    p = re.compile(regex_pattern)
    for m in p.finditer(seq_str):
        if m.start() == m.end():
            # even if there are groups in the lookahead, the first group should be the full match b/c that group surrounds the entire regex
            # so this will work whether or not there are groups in the lookahead
            match_seq = m.groups()[0]
        else:
            match_seq = seq_str[m.start() : m.end()]
        yield match_seq, m.start(), m.start() + len(match_seq)-1


