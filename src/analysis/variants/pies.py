'''
Creates the pie diagrams for the clinical data. No predictions are visualised here.
cwd is expected to be project root
'''

import pandas as pd
import os
import matplotlib.pyplot as plt
import re

# create a dataframe for the predictions
# it also contains coordinates to join with mmsplice output

df = pd.read_csv(os.path.join('variants', 'variants.csv'))

# determine SNV/MNV
df['Type'] = 'SNV'
df.loc[((df.REF.str.len() != 1) | (df.ALT.str.len() != 1)), 'Type'] = 'Not SNV'
df.loc[((df.REF.str.len() != 1) | (df.ALT.str.len() != 1)), 'Non-SNV Type'] = 'Substitution'
df.loc[(df.REF.str.len() > df.ALT.str.len()), 'Non-SNV Type'] = 'Deletion'
df.loc[(df.REF.str.len() < df.ALT.str.len()), 'Non-SNV Type'] = 'Insertion'

# Parse citations - remove citation chain from origin and create "multiple"
def removeCitationChain(source):
    return re.sub(r'\[[^\]]*\]', '', source).strip()

def clean_set(field):
    origin = set(field.split(', '))
    origin = set(map(lambda o: o.strip(), origin))
    return ', '.join(origin)

origin = df.Origin.apply(removeCitationChain)
origin = origin.apply(clean_set)
origin[origin.str.contains(',')] = 'Multiple'

# Count chromosomes
truncate_chroms_under = 80
chrom_count = ('chr'+df['#CHROM'].astype(str)).value_counts()
chrom_count_other = sum(chrom_count[chrom_count < truncate_chroms_under])
chrom_count = chrom_count[chrom_count >= truncate_chroms_under]
chrom_count['other'] = chrom_count_other

# Exact variant effect
exact_effect_subset = df.dropna(axis=0, subset=['AL', 'DL', 'AG', 'DG'], how='all')


# Pie plot
fig, axs = plt.subplots(2,4, figsize=(12,4))
axs = axs.reshape(-1)

standard_params = {
    'label': '',
}

origin.value_counts().plot(kind='pie', autopct='%1.0f%%', ax=axs[0], textprops={'fontsize': 8}, **standard_params)
df.SpliceAffecting.value_counts().plot(kind='pie', labels=['Affecting', 'Not Affecting'], autopct='%1.0f%%', ax=axs[1], startangle=70, **standard_params)
df.strand.value_counts().plot(kind='pie', autopct='%1.0f%%', ax=axs[2], startangle=70, **standard_params)
chrom_count.plot(kind='pie', autopct='%1.0f%%', ax=axs[3], **standard_params)
df.Type.value_counts().plot(kind='pie', autopct='%1.0f%%', explode=[0,.2], ax=axs[4], **standard_params)
df['Non-SNV Type'].value_counts().plot(kind='pie', autopct='%1.0f%%', ax=axs[5], **standard_params)
df.effect_type[df.effect_type != 'normal'].value_counts().plot(kind='pie', autopct='%1.0f%%', ax=axs[6], **standard_params)
exact_effect_subset.effect_type.value_counts().plot(kind='pie', autopct='%1.0f%%', ax=axs[7], **standard_params)


titles = ['A) Origin', 'B) Splice Affecting', 'C) Strand', 'D) Chromosome', 'E) Variant Type', 'F) Non-SNV Type', f'G) Broad Disruptions \n({len(df)} instances)', f'H) Exact Disruptions \n({len(exact_effect_subset)} instances)']

for ax, title in zip(axs, titles):
    ax.axis('equal')
    ax.set_title(title)

fig.tight_layout(h_pad=2)
fig.savefig(os.path.join('analysis', 'variants', 'pies.eps'), bbox_inches='tight')
fig.savefig(os.path.join('analysis', 'variants', 'pies.png'), bbox_inches='tight')