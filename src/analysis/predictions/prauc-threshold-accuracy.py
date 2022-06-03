'''
Creates the table holding PR-AUC, optimal threshold and acccuracy; creates the PR-AUC graph.
cwd is expected to be project root
'''

import pandas as pd
import os
import numpy as np
from sklearn.metrics import precision_recall_curve, auc, accuracy_score, roc_auc_score
import matplotlib.pyplot as plt
import aggregate_predictions

predictions, predictors = aggregate_predictions.read()

# ----- OPTIMAL THRESHOLDS -----

# calculate them
def calc_optimal_thresh(model_scores, labels):
    numerical = ~np.isnan(model_scores) & np.isfinite(model_scores)
    precision, recall, threshs = precision_recall_curve(labels[numerical], model_scores[numerical])

    accuracies_per_thresh = np.zeros_like(threshs)
    
    for i, thresh in enumerate(threshs):
        accuracies_per_thresh[i] = accuracy_score(labels, model_scores >= thresh)
    
    best_thresh = threshs[np.argmax(accuracies_per_thresh)]
    
    return best_thresh

   
optimal_thresholds = {}

for predictor in predictors:
    optimal_thresholds[predictor] = calc_optimal_thresh(predictions[predictor], predictions.SpliceAffecting)

print(f'Optimal thresholds: {optimal_thresholds}')

# ----- PR-AUC AND ACCURACY -----
def get_scores(data):
    table = pd.DataFrame()

    for predictor in predictors:
        coverage = (~pd.isna(data[predictor])).sum().astype(float) / len(data)

        # fill missing predictions with 0
        data_filled = data[predictor].fillna(0)
        data_filled[data_filled == np.inf] = 2 * max(data_filled[np.isfinite(data_filled)])
      
        accuracy_optimal = accuracy_score(data.SpliceAffecting, data_filled >= optimal_thresholds[predictor])

        precision, recall, _ = precision_recall_curve(data.SpliceAffecting, data_filled)
        pr_auc = auc(recall, precision)
        roc_auc = roc_auc_score(data.SpliceAffecting, data_filled)
        
        table = pd.concat([table, pd.Series({
            'Coverage': '%.0f%%' % (coverage*100),
            'AUC-PR': '%.2f%%' % (pr_auc*100),
            'AUC-ROC': '%.2f%%' % (roc_auc*100),
            'Optimal Threshold': '%.3f' % optimal_thresholds[predictor],
            'Accuracy': '%.2f%%' % (accuracy_optimal*100),
        }, name=predictor).to_frame().T])

        table.index.name = 'Algorithm'


    return table.sort_values(['AUC-PR', 'AUC-ROC', 'Accuracy'], ascending=True)

table = get_scores(predictions).filter(['Coverage', 'AUC-PR', 'AUC-ROC', 'Optimal Threshold', 'Accuracy'])
print(table)
# print(table.to_latex())
table.to_latex(os.path.join('analysis', 'predictions', 'prauc-threshold-accuracy.tex'))
table.to_csv(os.path.join('analysis', 'predictions', 'prauc-threshold-accuracy.csv'))
table.to_markdown(os.path.join('analysis', 'predictions', 'prauc-threshold-accuracy.md'))


# Calc PR-Curves
def p_r_auc(predictor, data):
    # fill nan with 0 and inf with 2*max
    data_filled = data[predictor].fillna(0)
    data_filled[data_filled == np.inf] = 2 * max(data_filled[np.isfinite(data_filled)])
    
    precision, recall, _ = precision_recall_curve(data.SpliceAffecting, data_filled)

    return {
        'precision': precision,
        'recall': recall,
        'auc': auc(recall, precision)
    }

p_r_aucs = pd.DataFrame({
    predictor: p_r_auc(predictor, predictions)
    for predictor in predictors
})

# Visualise them

# hack a second subplot to hold the figlegend in; this can be done better probably
fig, axs = plt.subplots(2, figsize=(4,5.4), gridspec_kw={'height_ratios': [9, 4]})
ax = axs[0]

# sort by AUC to make it easier to read
for _, row in p_r_aucs.T.sort_values('auc').iterrows():
    ax.step(row.recall, row.precision, label=f'{row.name} ({row.auc:.2f})')

# ax.axis('equal')
ax.set_xlim(-.02, 1.02)
ax.set_ylim(.5, 1.02)
ax.set_xlabel('Recall')
ax.set_ylabel('Precision')

ax = axs[1]
ax.axis('off')
fig.legend(loc='lower center')

fig.tight_layout()
fig.savefig(os.path.join('analysis', 'predictions', 'pr-auc.eps'))
fig.savefig(os.path.join('analysis', 'predictions', 'pr-auc.png'))