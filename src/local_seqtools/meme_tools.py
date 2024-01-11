from Bio import AlignIO, Seq, SeqIO
from Bio.motifs import meme


def import_meme_results_human_xml(filename):
    '''
    imports meme results from xml file
    returns list of motifs that have an instance in a human sequence.
    returned motifs are meme.Motif objects
    Checks for duplicate motifs and motifs occuring multiple times in the same sequence.
    this shouldn't happen unless meme was run with mod=anr
    I think I should rewrite this function if I want to parse a meme file with mod=anr

    '''
    with open(filename) as f:
        record = meme.read(f)
    mots_list = []
    dup_names = []
    mot_names = []
    for motif in record:
        dup_check = set()
        for inst in motif.instances:
            if inst.sequence_name in dup_check:
                print('DUPLICATE')
                print(inst.sequence_name)
                raise ValueError('duplicate')
            dup_check.add(inst.sequence_name)
            # print(inst.start, inst.motif_name, inst)
            if 'Homo_sapiens' in inst.sequence_name:
                # if inst.motif_name in mot_names:
                    # dup_names.append(inst.motif_name)
                # mot_names.append(inst.motif_name)
                mots_list.append(motif)
    # assert len(dup_names) == 0, f'there are duplicate motifs: {filename}\n{mot_names}'
    return mots_list


def import_meme_results_anr_human_xml(filename):
    '''
    imports meme results from xml file
    returns list of motifs that have an instance in a human sequence.
    returned motifs are meme.Motif objects
    '''
    with open(filename) as f:
        record = meme.read(f)
    mots_list = []
    mot_names = []
    for motif in record:
        dup_check = set()
        for inst in motif.instances:
            if 'Homo_sapiens' in inst.sequence_name:
                mots_list.append(motif)
                break
    return mots_list

    
def meme_motif_2_aln(motif):
    mot_aln_list = [SeqIO.SeqRecord(Seq.Seq(str(i)), id = i.sequence_name) for i in motif.instances]
    return AlignIO.MultipleSeqAlignment(mot_aln_list)

def meme_motif_2_posdict(motif):
    return {i.sequence_name:[i.start-1, i.start + len(str(i))-1-1] for i in motif.instances}

def meme_motif_2_posdict_anr(motif):
    pos = {i.sequence_name:[] for i in motif.instances}
    for i in motif.instances:
        pos[i.sequence_name].append([i.start-1, i.start + len(str(i))-1-1])
    return pos
