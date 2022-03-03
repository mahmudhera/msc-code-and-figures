# first, format the guidescan file properly
from Bio import trie
import pandas as pd


def complement(seq):
    complement_char = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = [complement_char[base] for base in bases]
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])


def compare(guides_kr, guides_gs):
    """
    :param guides_kr:list of tuples: (23-nt,strand)
    :param guides_gs:list of tuples: (20-nt,strand)
    :return: a minus b, a intersect b, b minus a
    """
    trie_kr = trie.trie()
    for guide in guides_kr:
        if guide[1] == '+':
            trie_kr[guide[0][:-3]] = 1
        else:
            trie_kr[reverse_complement(guide[0])[:-3]] = 1

    trie_gs = trie.trie()
    for guide in guides_gs:
        if guide[1] == '+':
            trie_gs[guide[0]] = 1
        else:
            trie_gs[reverse_complement(guide[0])] = 1

    intersection = 0
    for guide in trie_kr.keys():
        try:
            if trie_gs[guide] == 1:
                intersection += 1
        except:
            continue

    return len(trie_kr.keys())-intersection, intersection, len(trie_gs.keys())-intersection

guide_seqs_kr_p = pd.read_csv('krispmer-yeast.csv').tgt_in_plus.tolist()
guide_seqs_kr_n = pd.read_csv('krispmer-yeast.csv').tgt_in_minus.tolist()

guide_seqs_gs = pd.read_csv('GuideScan_batch_output.csv').gRNA.tolist()

intersection = 0
for guide in guide_seqs_gs:
    for longer_guide in guide_seqs_kr_n+guide_seqs_kr_p:
        if longer_guide.find(guide) >= 0:
            intersection += 1
            break

int2 = 0
for longer_guide in guide_seqs_kr_p:
    for guide in guide_seqs_gs:
        if longer_guide.find(guide) >= 0:
            int2 += 1
            break
        if reverse_complement(longer_guide).find(guide) >= 0:
            int2 += 1
            break

min_ = min(intersection, int2)

print( len(guide_seqs_kr_p)-min_, min_, len(guide_seqs_gs)-min_ )