import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import random

standard_genetic_code = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'N': ['AAT', 'AAC'],
        'D': ['GAT', 'GAC'],
        'C': ['TGT', 'TGC'],
        'Q': ['CAA', 'CAG'],
        'E': ['GAA', 'GAG'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAT', 'CAC'],
        'I': ['ATT', 'ATC', 'ATA'],
        'L': ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
        'K': ['AAA', 'AAG'],
        'M': ['ATG'],
        'F': ['TTT', 'TTC'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'W': ['TGG'],
        'Y': ['TAT', 'TAC'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'Stop': ['TAA', 'TAG', 'TGA']
    }

# 当我们研究一组 codon pair 是否使用同一个密码子时，我们希望看到的是前一个密码子的 tRNA 是否能被下一个密码子继续使用
# 所以需要统计的是，下一个密码子是否会是它本身/isoaccepting codon 还是其他 codon
accepting_dict = {
    'TTA': ['TTA', 'TTG'],
    'CTT': ['CTT', 'CTC', 'CTA'],
    'CTC': ['CTC', 'CTT'],
    'CTA': ['CTA', 'CTG'],
    'ATT': ['ATT', 'ATC', 'ATA'],
    'ATC': ['ATC', 'ATT'],
    'GTT': ['GTT', 'GTC', 'GTA'],
    'GTC': ['GTC', 'GTT'],
    'GTA': ['GTA', 'GTG'],
    'TCT': ['TCT', 'TCC', 'TCA'],
    'TCC': ['TCC', 'TCT'],
    'TCA': ['TCA', 'TCG'],
    'CCT': ['CCT', 'CCC', 'CCA'],
    'CCC': ['CCT', 'CCC'],
    'CCA': ['CCA', 'CCG'],
    'ACT': ['ACT', 'ACC', 'ACA'],
    'ACC': ['ACT', 'ACC'],
    'ACA': ['ACA', 'ACG'],
    'GCT': ['GCT', 'GCC', 'GCA'],
    'GCC': ['GCT', 'GCC'],
    'GCA': ['GCA', 'GCG'],
    'TAC': ['TAC', 'TAT'],
    'CAC': ['CAC', 'CAT'],
    'CAA': ['CAA', 'CAG'],
    'AAC': ['AAC', 'AAT'],
    'AAA': ['AAA', 'AAG'],
    'GAC': ['GAC', 'GAT'],
    'GAA': ['GAA', 'GAG'],
    'TGC': ['TGC', 'TGT'],
    'CGT': ['CGT', 'CGC', 'CGA'],
    'CGC': ['CGT', 'CGC'],
    'CGA': ['CGA', 'CGG'],
    'AGC': ['AGC', 'AGT'],
    'AGA': ['AGA', 'AGG'],
    'GGT': ['GGT', 'GGC', 'GGA'],
    'GGC': ['GGT', 'GGC'],
    'GGA': ['GGA', 'GGG'],
}


def codon_frequency(sequence):
    codon_count = {}
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in codon_count:
            codon_count[codon] += 1
        else:
            codon_count[codon] = 1
    return codon_count

def amino_acid_frequency(codon_count):
    amino_acid_count = {}
    for codon, count in codon_count.items():
        for amino_acid, codons in standard_genetic_code.items():
            if codon in codons:
                if amino_acid in amino_acid_count:
                    amino_acid_count[amino_acid] += count
                else:
                    amino_acid_count[amino_acid] = count
                break
    return amino_acid_count

def translate_sequence(sequence):
    amino_acid_sequence = ""
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        for amino_acid, codons in standard_genetic_code.items():
            if codon in codons:
                amino_acid_sequence += amino_acid
                break
    return amino_acid_sequence

def generate_random_sequence(amino_acid_seq, frequency):
    random_sequence = ""
    for amino_acid in amino_acid_seq:
        codons = standard_genetic_code[amino_acid]
        codon_weights = [frequency[codon] if codon in frequency else 0 for codon in codons]
        total_weight = sum(codon_weights)
        if total_weight == 0:
            chosen_codon = random.choice(codons)
        else:
            chosen_codon = random.choices(codons, weights=codon_weights, k=1)[0]
        random_sequence += chosen_codon
    return random_sequence

def generate_multiple_random_sequences(sequence, num_sequences=10):
    frequency = codon_frequency(sequence)
    amino_acid_seq = translate_sequence(sequence)
    random_sequences = [generate_random_sequence(amino_acid_seq, frequency) for _ in range(num_sequences)]
    return random_sequences

# 示例用法
sequence = "ATGCGACTACGATCGAGGGCCAT"
random_sequences = generate_multiple_random_sequences(sequence, 10)
for i, seq in enumerate(random_sequences):
    print(f"Random Sequence {i+1}: {seq}")