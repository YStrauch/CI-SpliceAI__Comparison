#!/bin/bash
source src/helpers.sh

# visualise clinical data, but no predictions
python src/analysis/variants/pies.py
python src/analysis/variants/distance-label.py

# visualise splice site predictions on CFTR (not variants)
python src/analysis/splicing/CFTR.py

# aggregate predictions
python src/analysis/predictions/aggregate_predictions.py

# and create all tables and figures
python src/analysis/predictions/prauc-threshold-accuracy.py
python src/analysis/predictions/predictive-errors.py
python src/analysis/predictions/variant-effects.py
python src/analysis/predictions/cis-fp-fn.py