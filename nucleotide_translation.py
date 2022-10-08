import gzip
import sys
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
                seq = line.rstrip()
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

genome_filename = sys.argv[1]
for header, seq in read_fasta(genome_filename):
    seq_copy = seq.upper()
    for shift in range(3):
        for start,stop in orf(seq_copy, shift):
            translate_preparation = [seq_copy[i:i+3] for i in range(start, stop, 3)]
            translate = map(lambda inf: codon_dictionary.get(inf, '*'), translate_preparation)
            dna = ''.join(translate)
            print(dna, '+')
    for shift in range(3):
        for start,stop in orf(swap_sequence(seq_copy), shift):
            translate_preparation = [seq_copy[i:i+3] for i in range(len(seq_copy) - 1 - stop, len(seq_copy) - 1 - start, 3)]
            translate = map(lambda inf: codon_dictionary.get(inf, '*'), translate_preparation)
            dna = ''.join(translate)
            print(dna, '-')
