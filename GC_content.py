import csv
import gzip
import matplotlib.pyplot as plt

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

genome_filename = "E_coli.fna"
for header, seq in read_fasta(genome_filename):
    GC_content = []
    with open('Proteins.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            Start, Stop, Strand = int(row['Start']), int(row['Stop']), row['Strand']
            if Strand == "+":
                seq_copy = seq[Start:Stop]
            else:
                seq_copy = swap_sequence(seq[Start:Stop])
            GC_content.append((seq_copy.count("G") + seq_copy.count("C")) / (Stop - Start + 1) * 100)
    plt.hist(GC_content, bins = 1000)
    plt.show()