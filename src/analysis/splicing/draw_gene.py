import numpy as np
import matplotlib.pyplot as plt
import re

def reduce_sequence(Y, threshold, num_replace = None):
    '''
    Removes recurring Y values when they occur `threshold` times, replaces with `num_replace` (defaults to `threshold`).
    Returns a tuple (Y_reduced, annotations) where annotations is a list of tuples (pos, num_removed)
    '''
    if num_replace is None:
        num_replace = threshold
    removed = []

    Y = ''.join(Y.astype(str))
    Y_processed = ''
    
    p = re.compile(r'(0{%d,})' % threshold)
    groups = p.split(Y)

    for i, group in enumerate(groups):
        if i % 2 == 0:
            Y_processed += group
        else:
            start = len(Y_processed)
            Y_processed += '0' * num_replace
            removed.append((start, len(group) - num_replace))

    return np.array(list(Y_processed), dtype=np.int8), removed

def apply_reduction_annotations(Y, annotations):
    for start, removed in annotations:
        Y = np.concatenate([Y[:start], Y[start+removed:]])
    return Y


def rolling_window(a, window):
    '''Rolls a window over an array a and emits a tuple of length len(window)'''
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def visualise_Y(ax: plt.Axes, Y: np.ndarray, plus_strand: bool):
    '''
        Visualises an array `Y` where:
            - 1 represents Acceptor
            - 2 represents Donor
            - 0 represents neither
    '''

    ax.axes.get_yaxis().set_visible(False)
    for spine in ['top', 'right', 'bottom', 'left']:
        ax.spines[spine].set_visible(False)

    sites = np.flatnonzero(Y)

    for site_index, (pos_before, pos, pos_after) in enumerate(rolling_window(np.pad(sites, 1, 'constant', constant_values=-1), 3)):
        site_before = Y[pos_before] if pos_before != -1 else None
        site = Y[pos]
        site_after = Y[pos_after] if pos_after != -1 else None

        starts_exon = site == 1
        starts_intron = site == 2
        cryptic_exon = site == site_after
        cryptic_intron = site == site_after

        # If an exon ends, draw the boundary *after* the nucleotide because all annotated nucs are within the exon
        if site == 2:
            pos += 1
        if Y[pos_after] == 2:
            pos_after += 1
        if Y[pos_before] == 2:
            pos_before += 1
        
        # 1. Draw vertical exon start/end
        style = 'solid'

        # Draw the vertical exon line solid only if it is not cryptic, i.e. 212 and 121
        if site_before is not None and site_before == site or site_after is not None and site_after == site:
            style='dotted'
        
        ax.vlines([pos], -1, 1, color='black', linestyles=style)


        # Exon and intron vertical lines are always drawn from left to right

        if starts_exon or cryptic_exon:
            style = 'solid' if starts_exon and not cryptic_exon else 'dotted'

            if site_after is None:
                # exon opens and is the last site, draw it until tx_end
                exon_width = len(Y) - pos
            else:
                exon_width = pos_after - pos

            for y in [-1,1]:
                ax.plot([pos, pos+(exon_width) ], [ y, y ], linestyle=style, color='black')
        
        if (starts_intron or cryptic_intron) and not site_index == len(sites) - 1:
            style = 'solid' if starts_intron and not cryptic_intron else 'dotted'

            if site_after is None:
                # draw intron until tx_end
                intron_width = len(Y) - pos
            else:
                intron_width = pos_after - pos

            ax.plot([pos, pos+(intron_width) ], [ 0, 0 ], linestyle=style, color='black')

        # Edge case: We start with an acceptor, need to draw to the left
        if site_index == 0 and starts_intron:
            for y in [-1,1]:
                ax.plot([0, pos], [ y, y ], linestyle='solid', color='black')

            
    # # Annotate removed parts
    # if truncate_annotations is not None:
    #     for pos, remove in truncate_annotations:
    #         ax.text(pos, 0, '-%d' % remove, fontsize=12, horizontalalignment='center', verticalalignment='bottom')

    # When setting xlim, the very left and right is truncated in the PDF, so add a small buffer
    buffer = max(1, int(len(Y) * .005))
    ax.set_xlim([-buffer, len(Y)+buffer])
    ax.set_ylim([-1.5, 1.5])

    ax.set_xticks([])


    if not plus_strand:
        ax.invert_xaxis()
    return ax


if __name__ == '__main__':
    Y = np.array(list('0000010010001'), dtype=np.int8)
    Y_reduced, annotations = reduce_sequence(Y, 3, 1)
    assert np.all(Y_reduced == np.array(list('0100101'), dtype=np.int8))
    assert annotations == [
        (0, 4),
        (5, 2)
    ]

    Y_reduced2 = apply_reduction_annotations(Y, annotations)
    assert (Y_reduced2 == Y_reduced).all()

