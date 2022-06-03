from typing import OrderedDict
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from pkg_resources import resource_filename

class Binning():
    def __init__(self):
        # we are using the splice sites provided by CI-SpliceAI
        self.annotation_table = pd.read_csv(resource_filename('cispliceai', os.path.join('data', 'grch38.csv'))).set_index('gene_id')

        self.acceptors = pd.concat([
            self.annotation_table.loc[self.annotation_table.strand == '+', 'jn_end'],
            self.annotation_table.loc[self.annotation_table.strand == '-', 'jn_start']
        ]).str.split(',').map(lambda lst: np.array(lst, dtype=int))

        self.donors = pd.concat([
            self.annotation_table.loc[self.annotation_table.strand == '+', 'jn_start'],
            self.annotation_table.loc[self.annotation_table.strand == '-', 'jn_end']
        ]).str.split(',').map(lambda lst: np.array(lst, dtype=int))

        self.df = pd.read_csv(os.path.join('variants', 'variants.csv')).set_index('ID')
        self.df.SpliceAffecting = self.df.SpliceAffecting.astype(bool)

        self.df = self.df.apply(self.annotate, axis=1)

        # split into acceptors and donors
        self.closest_is_acceptor = self.df[self.df['Closest SS Acceptor']]
        self.closest_is_donor = self.df[~self.df['Closest SS Acceptor']]

        # custom bins - tuple of (min,max), both inclusive
        bins = [
            (1,1),
            (2,2),
            (3,3),
            (4,4),
            (5,5),
            (6,10),
            (11,50),
            (51,200),
            (201,500),
            (501, np.inf)
        ]
        # negative offsets
        bins += [(-b[1], -b[0]) for b in bins]
        # and 0,0
        bins += [(0,0)]
        self.bins = sorted(bins, key=lambda key: key[0])

        # map bins to offset on the diagram
        self.x = np.arange(len(bins), dtype=int)

    def closest_ss_for_subset(self, var_pos: int, ref_len:int, sites: pd.Series):
        '''Returns the difference to the closest splice site in this subset of sites'''
        # if REF is more than one nucleotide, we check all positions between POS and POS+len(REF); this is not very efficient but simple enough
        closest_offset = np.inf

        for pos in range(var_pos, var_pos+ref_len):
            closest_per_gene = sites.map(lambda s: (pos - s)[np.argmin(np.abs(s-pos))])
            closest = closest_per_gene.iloc[np.argmin(closest_per_gene.abs())]
            if abs(closest) < abs(closest_offset):
                closest_offset = closest
        
        return closest_offset

    def closest_ss(self, pos: int, chrom:str, strand:str, ref_len: int):
        '''
            Returns the distance to the closest splice site and if it's an acceptor.
            +1 means closest SS is 1 downstream (downstream being relative to gene strand)
        '''
        # only look at genes on the correct chromosome and strand
        genes = self.annotation_table[
            (self.annotation_table.chr == chrom) &
            (self.annotation_table.strand == strand)
        ]

        closest_acceptor, closest_donor = self.closest_ss_for_subset(pos, ref_len, self.acceptors.loc[genes.index]), self.closest_ss_for_subset(pos, ref_len, self.donors.loc[genes.index])
        if strand == '-':
            closest_acceptor, closest_donor = -closest_acceptor, -closest_donor

        closest_is_acceptor = abs(closest_acceptor) < abs(closest_donor)
        closest_site = closest_acceptor if closest_is_acceptor else closest_donor

        return closest_site, closest_is_acceptor

    def annotate(self, row):
        row['Closest SS Offset'], row['Closest SS Acceptor'] = self.closest_ss(row.POS, 'chr'+row['#CHROM'], row.strand, len(row.REF))
        return row

    def signed(self, number: int):
        if number == 0:
            return number
        return f'{"+" if number >= 0 else "-"}{abs(number)}'

    def bin_label(self, bin):
        # create bin label
        label = self.signed(bin[0])
        if bin[0] != bin[1]:
            if np.isposinf(bin[1]):
                label = '$\geq$' + label
            elif np.isneginf(bin[0]):
                label = '$\leq$' + self.signed(bin[1])
            else:
                label += f' to {self.signed(bin[1])}'
        
        return label

    def bin(self, series, bins):
        counts = OrderedDict()
        for b in bins:
            data = series[(series >= b[0]) & (series <= b[1])]
            # if data.empty:
            #     continue

            counts[b] = len(data)
        return counts


    def draw_intron(self, ax: plt.Axes, start, sep_pos, stop, draw_y_start, draw_height, sep_offset=.1, sep_len=.1):
        # first intron line
        ax.hlines([draw_y_start+draw_height*.5], [start], [sep_pos-sep_len], color='black')

        # separator
        ax.plot([sep_pos-sep_len-sep_offset, sep_pos+sep_len-sep_offset], [draw_y_start+draw_height*.25, draw_y_start+draw_height*.75], color='black')
        ax.plot([sep_pos-sep_len+sep_offset, sep_pos+sep_len+sep_offset], [draw_y_start+draw_height*.25, draw_y_start+draw_height*.75], color='black')

        # second intron line
        ax.hlines([draw_y_start+draw_height*.5], [sep_pos+sep_len], [stop], color='black')

        # vertical boundary lines
        ax.vlines([start, stop], [draw_y_start, draw_y_start], [draw_y_start+draw_height, draw_y_start+draw_height], color='black')


    def draw_exon(self, ax: plt.Axes, start, stop, draw_y_start, draw_height):
        # two horizontal exon lines
        ax.hlines([draw_y_start, draw_y_start+draw_height], [start], [stop], color='black')



    # draw bin counts
    def draw_bin_counts(self, ax, series, offset, label=None, **barargs):
        counts = self.bin(series, self.bins)
        # x, y = np.array(list(range(offset, offset+len(counts)))), np.array(list(counts.values()))
        x = self.x + offset
        y = np.array(list(counts.values()))
        ax.bar(x, y, label=label, **barargs)
        return counts