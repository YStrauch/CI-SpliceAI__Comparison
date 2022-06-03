'''
Compares predictions of SpliceAI and CI-SpliceAI (trained on TRAIN models).
The CI-SpliceAI instance was trained on the TRAIN fold and copied into here.
Therefore both algorithms did not see this gene before.
'''

from pkg_resources import resource_filename
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ml_data
import draw_gene
import keras
from cispliceai.model import CISpliceAI

annotation_table = pd.read_csv(resource_filename('cispliceai', os.path.join('data', 'grch38.csv'))).set_index('gene_id')

annotation_table.jn_start = annotation_table.jn_start.map(lambda i: list(map(int, i.split(','))))
annotation_table.jn_end = annotation_table.jn_end.map(lambda i: list(map(int, i.split(','))))
gene = annotation_table.loc['ENSG00000001626'] # CFTR gene

cis_all = CISpliceAI()
cis_train = CISpliceAI(os.path.join('src', 'CI-SpliceAI-TRAIN.pb'))

def pred_from_h5(model_paths, X):
    models = [keras.models.load_model(path) for path in model_paths]
    preds = [model.predict(np.array([X]))[0] for model in models]
    # ensemble average
    preds = np.mean(preds, axis=0)
    return np.argmax(preds, axis=1)

def pred_CISpliceAI(model, X):
    preds = model.predict([X])[0]
    return np.argmax(preds, axis=1)

# Create ground truth & precictions

X, Y = ml_data.create_x(gene), ml_data.create_y(gene)
# add 5k of context to both sides
X = np.concatenate([np.zeros([5000,4]), X, np.zeros([5000,4])])

# predict them
Y_sai = pred_from_h5([resource_filename('spliceai', os.path.join('models', 'spliceai%d.h5')) % i for i in range(1,6)], X)
Y_cis_train= pred_CISpliceAI(cis_train, X)
Y_cis_all= pred_CISpliceAI(cis_all, X)

# Add an transcript start/end so the boxes are neatly closed
for y in [Y, Y_sai, Y_cis_train, Y_cis_all]:
    y[0] = 1
    y[-1] = 2

# truncate long regions where nothing happens so we can "zoom in" and see details
# first, combine all annotations
Y_shared = Y + Y_sai + Y_cis_train + Y_cis_all
# then, truncate sequences all in the same way
Y_shared_reduced, annotations = draw_gene.reduce_sequence(Y_shared, 200)
Y = draw_gene.apply_reduction_annotations(Y, annotations)
Y_sai = draw_gene.apply_reduction_annotations(Y_sai, annotations)
Y_cis_train = draw_gene.apply_reduction_annotations(Y_cis_train, annotations)
Y_cis_all = draw_gene.apply_reduction_annotations(Y_cis_all, annotations)

# Plot
fig, axs = plt.subplots(4, 1, figsize=(15, 4), sharex=True)

draw_gene.visualise_Y(axs[0], Y, gene.strand == '+')
draw_gene.visualise_Y(axs[1], Y_sai, gene.strand == '+')
draw_gene.visualise_Y(axs[2], Y_cis_train, gene.strand == '+')
draw_gene.visualise_Y(axs[3], Y_cis_all, gene.strand == '+')

# Mark incorrect predictions
for ax, Y_pred in zip(axs[1:], [Y_sai, Y_cis_train, Y_cis_all]):
    # add a red X to signal incorrect site
    for site in np.flatnonzero(Y + Y_pred): # iterate over all sites in this prediction and the ground truth
        if Y[site] != Y_pred[site]:
            # hacking y offset with -.1 so X aligns better
            ax.text(site, -.1, 'X', color='red', ha='center', va='center', weight='black', fontsize=14)

# Add tick labels
axs[-1].set_xticks(np.flatnonzero(Y_shared_reduced))
labels = np.flatnonzero(Y_shared) + gene.start
axs[-1].set_xticklabels(labels, rotation=90)


axs[0].set_title('GENCODE v37 GRCh38 (Collapsed)')
axs[1].set_title('SpliceAI')
axs[2].set_title('CI-SpliceAI (train chromosomes)')
axs[3].set_title('CI-SpliceAI')

for ax in axs:
    ax.set_ylim(-1.1, 1.1) # prevent lines being cropped
    ax.tick_params(length=0) # hide ticks themselves

fig.tight_layout()
fig.savefig(os.path.join('analysis', 'splicing', 'CFTR.eps'))
fig.savefig(os.path.join('analysis', 'splicing', 'CFTR.png'))

# Output stats
with open(os.path.join('analysis', 'splicing', 'stats.txt'), 'w') as f:
    f.write(f'In total, SpliceAI mispredicts {(Y != Y_sai).sum()} sites; CI-SpliceAI (ALL) mispredicts {(Y != Y_cis_all).sum()} sites; CI-SpliceAI (TRAIN) mispredicts {(Y != Y_cis_train).sum()} sites\n')
    f.write(f'SpliceAI has {np.logical_and(Y != 0, Y_sai == 0).sum()} false negatives and {np.logical_and(Y == 0, Y_sai != 0).sum()} false positives\n')
    f.write(f'CI-SpliceAI (ALL) has {np.logical_and(Y != 0, Y_cis_all == 0).sum()} false negatives and {np.logical_and(Y == 0, Y_cis_all != 0).sum()} false positives\n')
    f.write(f'CI-SpliceAI (TRAIN) has {np.logical_and(Y != 0, Y_cis_train == 0).sum()} false negatives and {np.logical_and(Y == 0, Y_cis_train != 0).sum()} false positives\n')
