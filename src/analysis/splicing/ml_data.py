import os
import numpy as np
import pyfaidx

# one-hot configurations
OH_X = np.asarray([[0, 0, 0, 0], # N
                   [1, 0, 0, 0], # A
                   [0, 1, 0, 0], # C
                   [0, 0, 1, 0], # G
                   [0, 0, 0, 1]],# T
                   dtype=bool)

# OH_Y = {
#     'neither':  np.asarray([1, 0, 0], dtype=bool),
#     'acceptor': np.asarray([0, 1, 0], dtype=bool),
#     'donor':    np.asarray([0, 0, 1], dtype=bool),
# }

def create_x(gene):
    '''Extracts sequence for a gene and encodes it one-hot. Reverse-complements if needed.'''
    fasta = pyfaidx.Fasta(os.path.join('third-party', 'hg', 'hg38.fa'))
    # -1 for start because fasta uses 0-index; not -1 for end to include the last nucleotide of the gene
    seq = fasta[gene.chr][gene.start-1:gene.end].seq

    # to numeric - reverse-complement if needed
    if gene.strand == '+':
        num = seq.upper().replace('N', '0').replace('A', '1').replace('C', '2').replace('G', '3').replace('T', '4')
    else:
        num = seq.upper().replace('A', '4').replace('C', '3').replace('G', '2').replace('T', '1').replace('N', '0')[::-1]

    num = np.asarray(list(map(int, list(num))))
    
    # one-hot encode
    oh = OH_X[num]
    return oh

def create_y(gene):
    '''Creates a one-hot encoded array with annotated ground truth (neither/acceptor/donor). Reverses if needed.'''
    if gene.strand == '+':
        donors, acceptors = gene.jn_start, gene.jn_end
    else:
        donors, acceptors = gene.jn_end, gene.jn_start

    y = np.array([0] * (gene.end - gene.start + 1)) # +1 to include last nucleotide of gene
    for site in donors:
        y[site - gene.start] = 2
    for site in acceptors:
        y[site - gene.start] = 1
    
    if gene.strand == '-':
        return y[::-1]
    return y