'''
Creates the variant effect comparison between all algorithms that support base-pair resolution of variant effects (i.e. SpliceAI, CI-SpliceAI, and VEP Sliding)
cwd is expected to be project root
'''

import pandas as pd
import aggregate_predictions
import os

predictions, predictors = aggregate_predictions.read()

def position_accuracy(kind, algo, data):
    # only look at data with annotated ground truth
    labelled = data[~data[kind].isna()].index
    pred = data.loc[labelled, f'{algo}_{kind}'] + predictions.POS[labelled]
    
    # calculate accuracy
    match = pred == data.loc[labelled, kind]
    accuracy = 100.0 * sum(match) / len(match)

    return f'{accuracy:.2f}%'

def create_table(algos):
    events = {
        'AG': 'Acceptor Gain',
        'AL': 'Acceptor Loss',
        'DG': 'Donor Gain',
        'DL': 'Donor Loss',
    }
    table = []

    for algo in algos:
        if f'{algo}_AL' not in predictions.columns:
            continue
        table.append(pd.Series({
            event: position_accuracy(event, algo, predictions)
            for event in events.keys()
        }, name=algo))
    
    table = pd.DataFrame(table).sort_values(list(events.keys()))[events.keys()]

    # table = table.rename(graph_labels, axis=0)
    table = table.rename(events, axis=1)
    table.index.name = 'Algorithm'
    
    return table
        

table = create_table(predictors)
print(table)
table.to_latex(os.path.join('analysis', 'predictions', 'variant-effects.tex'))
table.to_markdown(os.path.join('analysis', 'predictions', 'variant-effects.md'))