from typing import Tuple
import pysam
import os
import numpy as np
import sys
from maxentpy import maxent

sys.path.append(os.path.abspath('.'))
from src import helpers


genome_path = os.path.join('third-party', 'hg', 'hg38.fa')
# MES-specific number of nucleotides
n_nucs = {
    'donor': 9,
    'acceptor': 23,
}

n_nucs_before_jn = {
    'donor': 3,
    'acceptor': 20,
}

# load MES
matrices = {
    'donor': maxent.load_matrix5(),
    'acceptor': maxent.load_matrix3(),
}
scorers = {
    'donor': maxent.score5,
    'acceptor': maxent.score3,
}


def get_seq(line: pysam.VariantRecord, context: Tuple[int], rc: bool, apply_var: bool):
    assert len(line.alts) == 1, 'This code can only handle one alt per variant'
    alt = line.alts[0]
    seq = helpers.extract_sequence(genome_path, line.chrom, line.pos-context[0], line.pos+context[1])
    assert seq[context[0]:context[0]+len(line.ref)] == line.ref.upper(), 'mismatching REF annotation'
    
    if apply_var:
        seq = seq[:context[0]]+alt+seq[context[0]+len(line.ref):]
    
    if rc:
        return helpers.reverse_complement(seq)
    return seq


def pred_MES(line: pysam.VariantRecord):
    # runs MES on both strands around the variant and annotates the delta score
    annotations = []

    pred_lens = np.array(list(n_nucs.values())) + max(len(line.ref), len(line.alts[0])) - 1

    for reverse_strand in (False, True):
        DS_AL = DS_AG = DP_AL = DP_DG = DS_DG = DP_DG = 0

        for ss_type_num, (ss_type, nuc_len) in enumerate(n_nucs.items()): # site_type is acceptor/donor, each have their own length
            matrix = matrices[ss_type]
            scorer = scorers[ss_type]
            context = (pred_lens[ss_type_num]-1, pred_lens[ss_type_num])

            # predict ref and var
            preds = []
            for apply_var in (False, True):
                seq = get_seq(line, context, reverse_strand, apply_var)

                preds.append([scorer(seq[i:i+nuc_len], matrix=matrix) for i in range(len(seq) - nuc_len + 1)])

            preds_ref, preds_var = preds
            preds_ref, preds_var = np.array(preds_ref), np.array(preds_var)

            diff = abs(len(preds_ref) - len(preds_var))
            if diff:
                # compensate if the variant has another length than the reference

                def shorten_predictions(preds, index, length):
                    assert length != 1, 'Length is 1'
                    assert index > 0, 'Index must be positive'
                    assert index+length < len(preds), 'Cannot shorten more than available'
                    return  np.concatenate([
                        preds[:index],
                        [np.max(preds[index:index+length], axis=0)],
                        preds[index+length:]
                    ])

                def pad_predictions(preds, index, length):
                    assert length > 0, 'Length must be > 0'
                    assert index > 0, 'Index must be positive'
                    return  np.concatenate([
                        preds[:index],
                        np.zeros((length)),
                        preds[index:]
                    ])


                if len(preds_ref) > len(preds_var):
                    # DELETION, i.e. GAAG -> G. diff would be 3
                    # Pad variant
                    preds_var = pad_predictions(preds_var, context[0], diff)
                else:
                    # INSERTION, i.e. G -> GT. diff would be 1
                    # Shorten variant
                    preds_var = shorten_predictions(preds_var, context[0], diff+1)

                # Make sure we correctly sliced this madness
                assert len(preds_ref) == len(preds_var), 'Indel compensation failed'
            
            if reverse_strand:
                preds_ref = preds_ref[::-1]
                preds_var = preds_var[::-1]
            
            delta = preds_var - preds_ref

            ds_gain = np.nanmax(delta)
            ds_loss = -np.nanmin(delta)
            dp_gain = np.nanargmax(delta) - context[0]
            dp_loss = np.nanargmin(delta) - context[0]

            if ss_type == 'donor':
                DS_DL = ds_loss
                DP_DL = dp_loss
                DS_DG = ds_gain
                DP_DG = dp_gain
            else:
                DS_AL = ds_loss
                DP_AL = dp_loss
                DS_AG = ds_gain
                DP_AG = dp_gain

        annotations.append(f'{line.alts[0]}|{"-" if reverse_strand else "+"}|{DS_AG:.2f}|{DS_AL:.2f}|{DS_DG:.2f}|{DS_DL:.2f}|{DP_AG}|{DP_AL}|{DP_DG}|{DP_DL}')
    return annotations

# load variants
with pysam.VariantFile(os.path.join('variants', 'variants.vcf'), 'r') as file_in:
    header = file_in.header
    vcf = [v for v in file_in]

header.add_line(
    '##INFO=<ID=MES_SLIDING,Number=.,Type=String,Description="MES_SLIDING annotations. ALLELE|STRAND|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL (DS=Delta Score, DP=Delta Position, AG/AL=Acceptor Gain/Acceptor Loss, DG/DL=Donor Gain/Donor Loss.'
)

for line in vcf:
    line.info['MES_SLIDING'] = pred_MES(line)

with pysam.VariantFile(os.path.join('predictions', 'mes_sliding.vcf'), 'w', header=header) as file_out:
    for line in vcf:
        file_out.write(line)
