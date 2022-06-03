# this script is executed by setup.sh. cwd is expected to be project root

import pandas as pd
import os

gencode_gtf = os.path.join('third-party', 'mmsplice', 'gencode.gtf')
processed_gencode_gtf = os.path.join('third-party', 'mmsplice', 'protein_coding_no_dups.gtf')

# Pre-process gencode.gtf as described in https://kipoi.org/models/MMSplice/
if not os.path.isfile(processed_gencode_gtf):
    gtf = pd.read_csv(gencode_gtf, sep='\t', header=None, skiprows=5)
    gtf.columns = ['chrom', 'source', 'type', 'start', 'stop', '?', 'strand', '?', 'info']

    # filter to chromosomes 1-22 and X,Y only
    allowed = pd.Series(['chr%s' % d for d in list(range(1, 23)) + ['X', 'Y']])
    gtf = gtf[gtf.chrom.isin(allowed)]

    # drop duplicate exons
    is_exon = gtf['type'] == 'exon'
    exons = gtf[is_exon].drop_duplicates(['chrom', 'type', 'start', 'stop', 'strand'])
    gtf = gtf[~is_exon].append(exons)

    # filter genes to proteincoding ones
    gtf = gtf[(gtf.type != 'gene') | ((gtf.type == 'gene') & (gtf['info'].str.contains('protein_coding')))]

    gtf.to_csv(processed_gencode_gtf, sep='\t', header=False, index=False)
    del gtf

    # clean up gencode.gtf
    os.unlink(gencode_gtf)

