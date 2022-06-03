'''
Creates the distance-label plot which uses custom binning logic to show distributions for annotated variants and if they affect splicing per-offset around each site
cwd is expected to be project root
'''

import os
import matplotlib.pyplot as plt
import sys

sys.path.append(os.path.abspath('.'))
from src.binning import Binning

bins = Binning()

fig, ax = plt.subplots(figsize=(15, 5))

# first draw donor, then acceptor
donor_counts = bins.draw_bin_counts(ax, bins.closest_is_donor['Closest SS Offset'], 0, 'Closest Site is Donor')
acceptor_offset = len(donor_counts)
acceptor_counts = bins.draw_bin_counts(ax, bins.closest_is_acceptor['Closest SS Offset'], acceptor_offset, 'Closest Site is Acceptor')

assert sum(donor_counts.values())+sum(acceptor_counts.values()) == len(bins.df), 'Not all data points landed in a bin!'

# draw intron/exon doodles
draw_height = ax.get_ylim()[1] * .075
draw_start = ax.get_ylim()[1] + draw_height
ax.set_ylim((ax.get_ylim()[0], draw_start+draw_height * 6))

# draw doodle
exon_ends = list(donor_counts.keys()).index((0,0))
exon_starts = list(acceptor_counts.keys()).index((0,0)) + acceptor_offset

bins.draw_exon(ax, -.5, exon_ends+.5, draw_start, draw_height)
bins.draw_intron(ax, exon_ends+.5, len(bins.x)-.5, exon_starts-.5, draw_start, draw_height)
bins.draw_exon(ax, exon_starts-.5, len(acceptor_counts)+len(donor_counts)-.5, draw_start, draw_height)

# draw splice affecting bars on top
donor_counts_affecting = bins.draw_bin_counts(ax, bins.closest_is_donor.loc[bins.closest_is_donor.SpliceAffecting, 'Closest SS Offset'], 0, 'Splice Affecting', color='black', hatch='xxxx', fill=None, lw=None, alpha=.5)
acceptor_counts_affecting = bins.draw_bin_counts(ax,bins. closest_is_acceptor.loc[bins.closest_is_acceptor.SpliceAffecting, 'Closest SS Offset'], acceptor_offset, color='black', hatch='xxxx', fill=None, lw=None, alpha=.5)

# write the number of variants and the number of splice affecting variants on top of every bin
for x, count_all, count_affecting in zip(bins.x, donor_counts.values(), donor_counts_affecting.values()):
    ax.annotate('$\\frac{%d}{%d}$' % (count_affecting, count_all), xy=(x, count_all), textcoords='data', ha='center', va='bottom', fontsize=14)

for x, count_all, count_affecting in zip(bins.x+acceptor_offset, acceptor_counts.values(), acceptor_counts_affecting.values()):
    ax.annotate('$\\frac{%d}{%d}$' % (count_affecting, count_all), xy=(x, count_all), textcoords='data', ha='center', va='bottom', fontsize=14)

ax.set_xticks(range(len(donor_counts)+len(acceptor_counts)))
ax.set_xticklabels([bins.bin_label(b) for b in list(donor_counts.keys())+list(acceptor_counts.keys())], rotation=90)
# despine
for spine in ['right', 'top', 'left', 'bottom']:
    ax.spines[spine].set_visible(False)

# remove y axis
ax.set_yticks([])
# hide tick lines but keep tick labels
ax.tick_params(axis=u'both', which=u'both',length=0)

ax.legend(loc='upper right')
fig.tight_layout()
fig.savefig(os.path.join('analysis', 'variants', 'distance-label.eps'))
fig.savefig(os.path.join('analysis', 'variants', 'distance-label.png'))

# distribution measures
with open(os.path.join('analysis', 'variants', 'consensus_regions.txt'), 'w') as f:
    f.write('%d%% of sites are closer to a donor than an acceptor site\n' % (100*len(bins.closest_is_donor['Closest SS Offset'])/len(bins.df)))
    consensus_donor = donor_counts[(1,1)]+donor_counts[(2,2)]
    f.write('%d%% of donor sites are in consensus motif\n' % (100*consensus_donor/len(bins.closest_is_donor['Closest SS Offset'])))
    consensus_acceptor = acceptor_counts[(-1,-1)]+acceptor_counts[(-2,-2)]
    f.write('%d%% of acceptor sites are in consensus motif\n' % (100*consensus_acceptor/len(bins.closest_is_acceptor['Closest SS Offset'])))
    f.write('%d%% of all sites are in consensus motif\n' % (100*(consensus_donor+consensus_acceptor)/len(bins.df)))
