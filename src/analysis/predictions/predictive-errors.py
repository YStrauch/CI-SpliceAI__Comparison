'''
Creates the predictive error comparison between SpliceAI and CI-SpliceAI.
cwd is expected to be project root
'''

import os
import sys
import matplotlib.pyplot as plt
import aggregate_predictions
from matplotlib import cm

sys.path.append(os.path.abspath('.'))
from src.binning import Binning

cmap = cm.get_cmap('tab10') # default colour map

predictions, predictors = aggregate_predictions.read()

# Calculate predictive error and compare them between SpliceAI and CI-SpliceAI
error_spliceai = abs(predictions['SpliceAI']-predictions.SpliceAffecting)
error_cispliceai = abs(predictions['CI-SpliceAI']-predictions.SpliceAffecting)

def plot_bins(ax):

    # ----- BIN PLOT -----

    # Plot predictive error relative to closest splice site
    improved = predictions[error_spliceai > error_cispliceai]
    decreased = predictions[error_cispliceai > error_spliceai]
    same = predictions[error_cispliceai == error_spliceai]

    with open(os.path.join('analysis', 'predictions', 'predictive-error-change.txt'), 'w') as f:
        f.write(f'Predictions changed in {100*(len(predictions)-len(same))/len(predictions):.0f}% out of which {100*len(improved)/(len(predictions)-len(same)):.0f}% improved')

    bins = Binning()

    # first draw donor, then acceptor
    acceptor_offset = len(bins.x)

    # draw fn/fp bars
    donor_improved = bins.draw_bin_counts(ax, bins.closest_is_donor.filter(improved.index, axis=0)['Closest SS Offset'], -.2, 'Predictive Error Decreased', color=cmap.colors[0], width=.4)
    donor_decreased = bins.draw_bin_counts(ax, bins.closest_is_donor.filter(decreased.index, axis=0)['Closest SS Offset'], .2, 'Predictive Error Increased', color=cmap.colors[1], width=.4)
    acceptor_improved = bins.draw_bin_counts(ax, bins.closest_is_acceptor.filter(improved.index, axis=0)['Closest SS Offset'], acceptor_offset-.2, None, color=cmap.colors[0], width=.4)
    acceptor_decreased = bins.draw_bin_counts(ax, bins.closest_is_acceptor.filter(decreased.index, axis=0)['Closest SS Offset'], acceptor_offset+.2, None, color=cmap.colors[1], width=.4)

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
    for x, count_imp, count_dec in zip(bins.x, donor_improved.values(), donor_decreased.values()):
        if count_imp:
            ax.annotate(count_imp, xy=(x-.2, count_imp), textcoords='data', ha='center', va='bottom', fontsize=8)
        if count_dec:
            ax.annotate(count_dec, xy=(x+.2, count_dec), textcoords='data', ha='center', va='bottom', fontsize=8)

    for x, count_imp, count_dec in zip(bins.x, acceptor_improved.values(), acceptor_decreased.values()):
        if count_imp:
            ax.annotate(count_imp, xy=(x-.2+acceptor_offset, count_imp), textcoords='data', ha='center', va='bottom', fontsize=8)
        if count_dec:
            ax.annotate(count_dec, xy=(x+.2+acceptor_offset, count_dec), textcoords='data', ha='center', va='bottom', fontsize=8)

    ax.set_xticks(range(2*len(bins.x)))
    ax.set_xticklabels([bins.bin_label(b) for b in list(donor_improved.keys())+list(acceptor_improved.keys())], rotation=90)
    # despine
    for spine in ['right', 'top', 'left', 'bottom']:
        ax.spines[spine].set_visible(False)

    # remove y axis
    ax.set_yticks([])
    # hide tick lines but keep tick labels
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.legend(loc='upper right')

def plot_scatter(ax):

    # ----- SCATTER PLOT -----
    ax.scatter(error_cispliceai, error_spliceai, color='black', alpha=.1, label='Prediction')
    ax.plot([-1,2], [-1,2], label='Equal Error')
    ax.set_xlabel('Predictive Error CI-SpliceAI')
    ax.set_ylabel('Predictive Error SpliceAI')
    # fig.legend()
    ax.axis('equal')


fig, axs = plt.subplots(1, 2, figsize=(20,5), gridspec_kw={'width_ratios': [4, 0.905]})

plot_bins(axs[0])
plot_scatter(axs[1]) 
axs[0].set_title('A) Positions')
axs[1].set_title('B) Magnitude')

fig.tight_layout()

# for some reason this needs to be called after tight_layout
axs[1].set_xlim(-.02, 1.02)
axs[1].set_ylim(-.02, 1.02)

fig.savefig(os.path.join('analysis', 'predictions', 'predictive-errors.png'))