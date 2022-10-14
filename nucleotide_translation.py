import random
import gzip
import matplotlib.pyplot as plt
codon_dictionary = {'TGA': '*', 'GCG': 'A', 'CGA': 'R', 'ATA': 'I', 'AGA': 'R', 'TAA': '*', 'TTT': 'F', 'GAG': 'E', 'CTT': 'L', 'CGT': 'R', 'CTC': 'L', 'CTG': 'L', 'TGT': 'C', 'CCA': 'P', 'AAT': 'N', 'GTC': 'V', 'GAC': 'D', 'GAT': 'D', 'TAT': 'Y', 'AAA': 'K', 'GTA': 'V', 'TAG': '*', 'CGC': 'R', 'GCA': 'A', 'TCG': 'S', 'GCT': 'A', 'GCC': 'A', 'TGG': 'W', 'TTC': 'F', 'CCC': 'P', 'TTG': 'L', 'CGG': 'R', 'GGC': 'G', 'AGG': 'R', 'TCC': 'S', 'CCT': 'P', 'GGT': 'G', 'GGG': 'G', 'TCA': 'S', 'AGC': 'S', 'CAG': 'Q', 'CAC': 'H', 'ATC': 'I', 'GAA': 'E', 'GTG': 'V', 'CCG': 'P', 'CAT': 'H', 'AAG': 'K', 'ATG': 'M', 'AAC': 'N', 'TAC': 'Y', 'TGC': 'C', 'CTA': 'L', 'TCT': 'S', 'ATT': 'I', 'ACG': 'T', 'AGT': 'S', 'GTT': 'V', 'TTA': 'L', 'CAA': 'Q', 'GGA': 'G', 'ACC': 'T', 'ACA': 'T', 'ACT': 'T'}

def read_fasta(path):
    header = None
    seq = None
    if path.endswith('.gz'):
        open_funk = gzip.open
    else:
        open_funk = open
    with open_funk(path, 'rt') as file:
        for line in file:
            if line[0] == '>':
                if header is not None:
                    yield(header, seq)
                header = line[1:-1]
                seq = ''
            else:
                seq += line.rstrip()
        if header is not None:
            yield(header, seq)

def swap_sequence(seq):
    reversed_seq = seq[::-1]
    transition = str.maketrans({
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'
    })
    return reversed_seq.translate(transition)

def orf(seq, shift = 0):
    start = '_'
    for number in range(shift, len(seq) - 2, 3):
        chunk = seq[number:number+3]
        if chunk == 'ATG' and start == '_':
            start = number
        elif codon_dictionary.get(chunk) == '*' and start != '_':
            yield(start, number + 2)
            start = '_'

def generate_nucleotids(length):
    letters = ['A', 'T', 'C', 'G']
    genom = ''
    for necleotid in range(length):
        chosen = random.randint(0, 3)
        genom += letters[chosen]
    return genom


genome_filename = "E_coli.fna"
for header, seq in read_fasta(genome_filename):
    for both in 1, 2:
        if both == 1:
            seq_copy = seq.upper()
            name = 'E_coli'
        else:
            seq_copy = generate_nucleotids(len(seq))
            name = 'Generated'
        Length, Amount = [], []
        Amino = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*']
        Nucl, NuclAmou, AminoAmou = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []], [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []], [0] * len(Amino)
        sum = 0
        for shift in range(3):
            for start,stop in orf(seq_copy, shift):
                if Length.count(stop - start + 1) == 0:
                    Length.append(stop - start + 1)
                    Amount.append(1)
                    sum += 1
                else:
                    Amount[Length.index(stop - start + 1)] += 1
                    sum += 1
                translate_preparation = []
                for i in range(start, stop, 3):
                    translate_preparation.append(seq_copy[i:i+3])
                    Shortcut = Amino.index(codon_dictionary.get(seq_copy[i:i+3]))
                    AminoAmou[Shortcut] += 1
                    if Nucl[Shortcut].count(seq_copy[i:i+3]) == 0:
                        Nucl[Shortcut].append(seq_copy[i:i+3])
                        NuclAmou[Shortcut].append(1)
                    else:
                        NuclAmou[Shortcut][Nucl[Shortcut].index(seq_copy[i:i+3])] += 1
                translate = map(lambda inf: codon_dictionary.get(inf, '*'), translate_preparation)
                dna = ''.join(translate)
        seq_copy = swap_sequence(seq_copy)
        for shift in range(3):
            for start,stop in orf(seq_copy, shift):
                if Length.count(stop - start + 1) == 0:
                    Length.append(stop - start + 1)
                    Amount.append(1)
                    sum += 1
                else:
                    Amount[Length.index(stop - start + 1)] += 1
                    sum += 1
                translate_preparation = []
                for i in range(start, stop, 3):
                    translate_preparation.append(seq_copy[i:i+3])
                    Shortcut = Amino.index(codon_dictionary.get(seq_copy[i:i+3]))
                    AminoAmou[Shortcut] += 1
                    if Nucl[Shortcut].count(seq_copy[i:i+3]) == 0:
                        Nucl[Shortcut].append(seq_copy[i:i+3])
                        NuclAmou[Shortcut].append(1)
                    else:
                        NuclAmou[Shortcut][Nucl[Shortcut].index(seq_copy[i:i+3])] += 1
                translate = map(lambda inf: codon_dictionary.get(inf, '*'), translate_preparation)
                dna = ''.join(translate)
        plt.bar(Amino, AminoAmou)
        plt.title('Кол-во раз аминокислоты были задествованы в ' + str(name))
        plt.xlabel('Названия аминокислот')
        plt.ylabel('Их кол-во появлений')
        plt.show()
        for number in range(len(Amino)):
            plt.bar(Nucl[number], NuclAmou[number])
            plt.title('Кол-во ра изпользования нуклеоидов ' + Amino[number])
            plt.xlabel('Названия нуклеоидов')
            plt.ylabel('Их кол-во раз использований')
            plt.show()
