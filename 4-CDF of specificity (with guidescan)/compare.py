from get_cfd_score import get_score
from generate_adjacent_mers import generate_adjacent_mers
import dna_jellyfish as jellyfish
from Bio import trie
import pandas as pd
import subprocess
from matplotlib import pyplot as plt
import numpy as np

max_hd = 2


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
    c1 = genome_file[mer]
    c2 = target_file[mer]
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


def generate_guidescan_specificity(guide, genome_jf):
    """
    generates the specificity of a guide sequence (including NGG-5' PAM)
    :param guide: 2-tuple of guide sequence (23-nts), and its strand type (either '+' or '-')
    :param genome_jf: jellyfish file of entire reference (generated for 23-nts)
    :return: floating point specificity score
    """
    gRNA_sequence = guide[0]
    gRNA_strand = guide[1]
    trie = generate_adjacent_mers(gRNA_sequence, max_hd)
    inverse_specificity = 0.0
    for seq in trie.keys():
        genome_count = get_count_in_genome(genome_jf, seq)
        cfd_score = generate_cfd_score(gRNA_sequence, seq, gRNA_strand)
        inverse_specificity += 1.0 * cfd_score * genome_count
    try:
        return 1.0 / inverse_specificity
    except:
        return 0.0


def annotate_guides_with_specificity(guides, genome_jf):
    """
    takes a list of guides and annotates them with specificity score
    :param guides: list of guides, which are 2-tuple of grna(23-nts) and strand('+' or '-')
    :param genome_jf: jellyfish filename (generated for 23-mers)
    :return: list
    """
    lst = [(guide, generate_guidescan_specificity(guide, genome_jf)) for guide in guides]
    return lst


def get_specificity(guidescan_output_filename, genome_fasta_filename, genome_jf_filename):
    """
    generates specificity of all guides in the gs output file and generates the specificity
    :param guidescan_output_filename: guidescan output filename
    :param genome_fasta_filename: self explanatory
    :param genome_jf_filename: jellyfish binary filename generated for the same genome, and for 23-mers
    :return: list of annotated guides
    """
    guides = generate_grna_sequenes(guidescan_output_filename, genome_fasta_filename)
    strands = pd.read_csv(guidescan_output_filename).strand.tolist()
    print (len(guides))
    print (len(strands))
    lst = annotate_guides_with_specificity(zip(guides, strands), genome_jf_filename)
    return lst


def get_specificity_scores_in_list(guidescan_output_filename, genome_fasta_filename, genome_jf_filename):
    """
    :return: self explanatory
    """
    lst = get_specificity(guidescan_output_filename, genome_fasta_filename, genome_jf_filename)
    return [s for (g, s) in lst]


def generate_plottable_points(specificity_scores):
    """
    with a list of specificity scores in a list, generates a CDF
    :param specificity_scores: list of floating point values
    :return: pair of lists, x and y, plotting which will plot the CDF
    """
    vals = np.arange(0, 1.01, 0.025)
    s = np.array(specificity_scores)
    ys = []
    for v in vals:
        ys.append(np.sum(np.count_nonzero(s <= v)))
    yss = [y / len(specificity_scores) for y in ys]
    return vals, yss


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


def generate_inverted_specificity_from_genome(guides, genome_jf_filename, target_jf_filename, target_count):
    """
    generates the inverted specificity score used in krispmer, but using the genome, not expectations
    :param guides: list of guides, which are 2-tuple of grna(23-nts) and strand('+' or '-')
    :param genome_jf_filename: jellyfish filename (generated for 23-mers)
    :param target_region_filename: the target string as fasta
    :return: dictionary containing the scores. dic[gRNA sequence] -> the inv_spec score
    """
    qf = jellyfish.QueryMerFile(genome_jf_filename)
    tgt_qf = jellyfish.QueryMerFile(target_jf_filename)
    dic = {}
    global max_hd
    for candidate, strand_type in guides:
        trie = generate_adjacent_mers(candidate, max_hd)
        val1 = 0.0
        val2 = 0.0
        for mer in trie.keys():
            merDNA = jellyfish.MerDNA(mer)
            revcompMerDNA = jellyfish.MerDNA(reverse_complement(mer))
            cutting_probability = generate_cfd_score(candidate, mer, '+')
            val1 += max(qf[merDNA], qf[revcompMerDNA]) * cutting_probability
            val2 += max(tgt_qf[merDNA], tgt_qf[revcompMerDNA]) * cutting_probability
            #print (mer + ' ' + str(cutting_probability))
        try:
            dic[candidate] = 1.0 * val1 / (val2 * target_count)
        except:
            dic[candidate] = -1
    return dic


def generate_inverted_specificity_scores_from_genome():
    df = pd.read_csv('krispmer-yeast.csv')
    tgt_regions_plus = df.tgt_in_plus.tolist()
    tgt_regions_minus = df.tgt_in_minus.tolist()
    strands = df.strand.tolist()
    guides = []
    for tgt_plus, tgt_minus, strand in zip(tgt_regions_plus, tgt_regions_minus, strands):
        if strand == '+':
            gRNA_seq = tgt_plus
            guides.append((gRNA_seq, strand))
        else:
            gRNA_seq = tgt_minus
            guides.append((gRNA_seq, strand))
    # guides = [k for k in guides if k[0] == 'TACGTCTTGGATTTCTACGGAGG']
    # print(guides)
    dik = generate_inverted_specificity_from_genome(guides, 'yeast_genome.jf', 'yeast_target.jf', 2)
    print(dik.values())


if __name__ == '__main__':
    l = get_specificity_scores_in_list('gs.csv', 'genome.fasta', 'jf_human_genome')
    print(l)