from get_cfd_score import get_score
from generate_adjacent_mers import generate_adjacent_mers
import dna_jellyfish as jellyfish
from Bio import trie
import pandas as pd
import subprocess
from matplotlib import pyplot as plt
import numpy as np

max_hd = 3


def complement(seq):
    """
    generates the complement sequence, e.g.: ACGT-->TCGA
    :param seq: dna sequence
    :return: complement
    """
    complement_char = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(str(seq))
    bases = [complement_char[base] for base in bases]
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])


def generate_cfd_score(gRNA, seq, strand):
    """
    generates Doench et. al. CDF score
    :param gRNA: guide RNA (23 nts)
    :param seq: sequence where the guide RNA is intended to make a cut (23nts)
    :param strand: +:5-3 strand
    :return: a floating points cdf score, indicating how likely the cut is
    """
    if strand == '+':
        return get_score(gRNA, seq)
    else:
        return get_score(reverse_complement(gRNA), reverse_complement(seq))


def generate_grna_sequenes(guidescan_filename, genome_filename):
    """
    guidescan output file only has genomic coordinates, with only the guide sequence (excluding PAM). this method
    generates the gRNA sequence with the PAM (23-nts)
    :param guidescan_filename: guidescan output file
    :param genome_filename: fasta filename of genome
    :return: list containing the guide sequences, in the same order of the guidescan_filename
    """
    df = pd.read_csv(guidescan_filename)
    start_coords = df.start.tolist()
    end_coords = df.end.tolist()
    chromosem = df.chromosome.tolist()
    grna_seqs = []
    for c, s, e in zip(chromosem, start_coords, end_coords):
        cmd = 'samtools faidx ' + genome_filename + ' ' + c + ':' + str(s) + '-' + str(e)
        args = cmd.split(' ')
        out = subprocess.check_output(args)
        grna_seqs.append(out.split('\n')[1])
    return grna_seqs


def get_num_off_targets(genome_jf, target_jf, target_count, sequence):
    """
    generates # of off-target cuts made by a guide sequence (with PAM, 23 nts)
    :param genome_jf: jellyfish file of entire genome (generated for 23-mers)
    :param target_jf: jellyfish file of the target sequence (generated for 23-mers)
    :param target_count: # of times a target string appears in the genome
    :param sequence: the gRNA sequence including NGG PAM (23 nts)
    :return: count (integer)
    """
    genome_file = jellyfish.QueryMerFile(genome_jf)
    target_file = jellyfish.QueryMerFile(target_jf)
    mer = jellyfish.MerDNA(sequence)
    rev_mer = jellyfish.MerDNA(reverse_complement(sequence))
    c1 = max(genome_file[mer], genome_file[rev_mer])
    c2 = max(target_file[mer], target_file[rev_mer])
    ot_count = max(c1 - target_count * c2, 0)
    return ot_count


def get_count_in_genome(genome_jf, sequence):
    """
    retrieves the # of times a sequence is in the genome
    :param genome_jf: jellyfish file of entire reference (generated for 23-nts)
    :param sequence: the sequence that is to be counted
    :return: count (integer)
    """
    genome_file = jellyfish.QueryMerFile(genome_jf)
    return max(genome_file[jellyfish.MerDNA(sequence)], genome_file[jellyfish.MerDNA(reverse_complement(sequence))])


def read_target_region(filename):
    """
    reads a fasta file and generates the target region as a string
    :param filename: fasta file (target region)
    :return: string without the line containing '>'
    """
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content][1:]
    return ''.join(content)


if __name__ == '__main__':
    # get the guides in krispmer
    # filter those with score <= 3
    # get the count with ot > 0
    df = pd.read_csv('krispmer-yeast.csv')
    guides = df.tgt_in_plus.tolist()
    count = 0
    for guide in guides:
        if get_num_off_targets('yeast_genome.jf', 'yeast_target.jf', 2, guide) > 0:
            count += 1
    print(count)
    gs_guides = generate_grna_sequenes('GuideScan_batch_output.csv', 'yeast.fasta')
    count = 0
    for guide in gs_guides:
        if get_num_off_targets('yeast_genome.jf', 'yeast_target.jf', 2, guide) > 0:
            count += 1
    print(count)