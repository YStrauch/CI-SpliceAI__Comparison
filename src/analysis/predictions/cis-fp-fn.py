'''
Creates the distance-label plot which uses custom binning logic to show distributions for annotated variants and if they affect splicing per-offset around each site
cwd is expected to be project root
'''

import os
import matplotlib.pyplot as plt
import sys
import pandas as pd
from matplotlib import cm
import aggregate_predictions

cmap = cm.get_cmap('tab10') # default colour map


sys.path.append(os.path.abspath('.'))
from src.binning import Binning

df, _ = aggregate_predictions.read()

# use the optimal threshold derived by prauc-threshold-accuracy analysis
threshold = pd.read_csv(os.path.join('analysis', 'predictions', 'prauc-threshold-accuracy.csv')).set_index('Algorithm')['Optimal Threshold']['CI-SpliceAI']

fp = df[(df['CI-SpliceAI'] >= threshold) & ~df.SpliceAffecting]
fn = df[(df['CI-SpliceAI'] < threshold) & df.SpliceAffecting]

bins = Binning()

fig, ax = plt.subplots(figsize=(15, 5))

# first draw donor, then acceptor
acceptor_offset = len(bins.x)

# draw fn/fp bars
donor_fp = bins.draw_bin_counts(ax, bins.closest_is_donor.filter(fp.index, axis=0)['Closest SS Offset'], -.2, 'False Positive', color=cmap.colors[0], width=.4)
donor_fn = bins.draw_bin_counts(ax, bins.closest_is_donor.filter(fn.index, axis=0)['Closest SS Offset'], .2, 'False Negative', color=cmap.colors[1], width=.4)
acceptor_fp = bins.draw_bin_counts(ax, bins.closest_is_acceptor.filter(fp.index, axis=0)['Closest SS Offset'], acceptor_offset-.2, None, color=cmap.colors[0], width=.4)
acceptor_fn = bins.draw_bin_counts(ax, bins.closest_is_acceptor.filter(fn.index, axis=0)['Closest SS Offset'], acceptor_offset+.2, None, color=cmap.colors[1], width=.4)

# draw intron/exon doodles
draw_height = ax.get_ylim()[1] * .075
draw_start = ax.get_ylim()[1] + draw_height*.5
ax.set_ylim((ax.get_ylim()[0], draw_start+draw_height * 4))

# draw doodle
exon_ends = list(bins.bins).index((0,0))
exon_starts = list(bins.bins).index((0,0)) + acceptor_offset

bins.draw_exon(ax, -.5, exon_ends+.5, draw_start, draw_height)
bins.draw_intron(ax, exon_ends+.5, len(bins.x)-.5, exon_starts-.5, draw_start, draw_height)
bins.draw_exon(ax, exon_starts-.5, 2*len(bins.x)-.5, draw_start, draw_height)


# add numerical values
for x, count_fp, count_fn in zip(bins.x, donor_fp.values(), donor_fn.values()):
    if count_fp:
        ax.annotate(count_fp, xy=(x-.2, count_fp), textcoords='data', ha='center', va='bottom', fontsize=10)
    if count_fn:
        ax.annotate(count_fn, xy=(x+.2, count_fn), textcoords='data', ha='center', va='bottom', fontsize=10)

for x, count_fp, count_fn in zip(bins.x, acceptor_fp.values(), acceptor_fn.values()):
    if count_fp:
        ax.annotate(count_fp, xy=(x-.2+acceptor_offset, count_fp), textcoords='data', ha='center', va='bottom', fontsize=10)
    if count_fn:
        ax.annotate(count_fn, xy=(x+.2+acceptor_offset, count_fn), textcoords='data', ha='center', va='bottom', fontsize=10)

ax.set_xticks(range(2*len(bins.x)))
ax.set_xticklabels([bins.bin_label(b) for b in list(donor_fp.keys())+list(acceptor_fp.keys())], rotation=90)
# despine
for spine in ['right', 'top', 'left', 'bottom']:
    ax.spines[spine].set_visible(False)

# remove y axis
ax.set_yticks([])
# hide tick lines but keep tick labels
ax.tick_params(axis=u'both', which=u'both',length=0)

ax.legend(loc='upper right')
fig.tight_layout()
fig.savefig(os.path.join('analysis', 'predictions', 'cis-fp-fn.eps'))
fig.savefig(os.path.join('analysis', 'predictions', 'cis-fp-fn.png'))