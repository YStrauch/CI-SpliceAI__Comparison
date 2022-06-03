import pandas as pd
import sys
import os
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import helpers

df = pd.read_csv(os.path.join('..', 'variants', 'variants.csv')).set_index('ID')
genome_path = sys.argv[1]

# -- CHECK DATA --

# Make sure that the reference nucleotides match
df = df.reset_index()

if np.any(df['REF'] == '-') or np.any(df['ALT'] == '-') or pd.isna(df['REF']).any() or pd.isna(df['ALT']).any():
    raise Exception('Missing REF or ALT annotations. Make sure to specify neighbouring nucleotides.')
    
df['validate_ref'] = helpers.extract_sequences(genome_path, df['#CHROM'].astype(str), df['POS'], df['POS'] + df.REF.str.len())

df = df.set_index('ID')
mismatched = df.REF != df.validate_ref

# Assert that there are no mismatches
if sum(mismatched) != 0:
    print(df[mismatched])
    raise Exception("Mismatch between REF annotations and genome in %d cases" % sum(mismatched))


# Append 'chr'
df['#CHROM'] = 'chr'+df['#CHROM']

# Add missing cols
df['FILTER'] = '.'
df['INFO'] = '.'

# Export in vcf format
csv = df.reset_index().filter(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'FILTER', 'INFO']).to_csv(index=False, sep='\t', header=True)

with open('plain.vcf', 'w') as f:
    f.writelines('''##fileformat=VCFv4.2
##reference=GRCh38\n''')
    f.write(csv)