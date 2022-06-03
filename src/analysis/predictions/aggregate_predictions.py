'''
Aggregates all predictions with ground truth into predictions/joined.
Assumes cwd to be project root
'''
import pandas as pd
import pysam
import os
import numpy as np
import sys

def aggregate():
    predictors = set()

    sys.path.append(os.path.abspath('.'))
    from src import helpers 

    # read in variant locations and ground truth
    df = pd.read_csv(os.path.join('variants', 'variants.csv')).set_index('ID')

    # ----- PARSE VCF PREDICTION FILES -----

    # functions to extract binary prediction and corresponding effect position (if possible) from VCF info column
    # each emits a tuple of (score, (pos_ag, pos_al, pos_dg, pos_dl)); for algorithms where positions are unknown we emit (pos, None)

    def parse_cis_format_annotation(info):
        '''
        Helper function used for SpliceAI, CI-SpliceAI and MES(Sliding).
        Will look at all annotations and extract the highest delta score, and the respective positions of the highest AG,AL,DG,DL
        '''
        if len(info) < 1:
            return None, None

        annotations = np.zeros((len(info), 4, 2), dtype=np.float32) # shape: no. annotations x AG/DG/AL/DL x score/position
        for i, annotation in enumerate(info):
            annotation = annotation.split('|')
            annotations[i,:,0] = annotation[2:2+4]
            annotations[i,:,1] = annotation[6:6+4]

        DP_AG = annotations[:,0,1][np.argmax(annotations[:,0,0])].astype(int)
        DP_AL = annotations[:,1,1][np.argmax(annotations[:,1,0])].astype(int)
        DP_DG = annotations[:,2,1][np.argmax(annotations[:,2,0])].astype(int)
        DP_DL = annotations[:,3,1][np.argmax(annotations[:,3,0])].astype(int)

        return np.max(annotations[:,:,0]), (DP_AG,DP_AL,DP_DG,DP_DL)

    def parse_cis_info(info):
        return parse_cis_format_annotation(info['CI-SpliceAI'])

    def parse_spliceai_info(info):
        if 'SpliceAI' not in info:
            return None, None
        
        annotations = info['SpliceAI']
        # filter out annotations with '.' as score
        annotations = [annotation for annotation in annotations if '.' not in annotation.split('|')[2:6]]
        return parse_cis_format_annotation(annotations)

    def parse_mms_sliding_info(info):
        return parse_cis_format_annotation(info['MES_SLIDING'])


    def parse_mes_vep_info(info):
        max_annotation = -1
        for annotation in info['CSQ']:
            diff = annotation.split('|')[-2] # second last entry is "MaxEntScan_diff"

            if diff == '' or diff == '.':
                continue
            max_annotation = max(max_annotation, abs(float(diff)))
        
        return (None, None) if max_annotation == -1 else (max_annotation, None)

    def parse_squirls_info(info):
        annotation = info['SQUIRLS_SCORE'][0]
        max_annotation = -1

        for _, score in map(lambda field: field.split('='), annotation.split('|')[1:]):
            max_annotation = max(max_annotation, abs(float(score)))

        return (None, None) if max_annotation == -1 else (max_annotation, None)



    def append_from_vcf(predictions, fname, col_name, extract_fn):
        # predictions[col_name] = None
        with pysam.VariantFile(os.path.join('predictions', fname)) as f:
            for row in f:
                score, positions = extract_fn(row.info)
                predictions.loc[row.id,col_name] = score
                if positions is not None:
                    for position, identifier in zip(positions, ['AG','AL','DG','DL']):
                        df.loc[row.id,f'{col_name}_{identifier}'] = position

        predictors.add(col_name)
        return predictions

    df = append_from_vcf(df, 'cis.vcf', 'CI-SpliceAI', parse_cis_info)
    df = append_from_vcf(df, 'spliceai.vcf', 'SpliceAI', parse_spliceai_info)
    df = append_from_vcf(df, 'mes_vep.vcf', 'MES (VEP)', parse_mes_vep_info)
    df = append_from_vcf(df, 'mes_sliding.vcf', 'MES (Sliding)', parse_mms_sliding_info)
    df = append_from_vcf(df, 'squirls.vcf', 'SQUIRLS', parse_squirls_info)



    # ----- PARSE MMSplice TSV FILES -----
    # they don't write variant IDs, so we need to join data based on POS/REF/ALT

    def append_from_mms_tsv(predictions, predictor, col_name):
        f_in = os.path.join('predictions', f'mmsplice_{predictor}.tsv')

        # maps tsv values from output tsv to the col names from our csv
        colmap = {
            'metadata/variant/CHROM': '#CHROM',
            'metadata/variant/POS': 'POS',
            'metadata/variant/REF': 'REF',
            'metadata/variant/ALT': 'ALT',
            'metadata/variant/POS': 'POS',
            'preds': col_name
        }


        index_cols = [col for col in colmap.values() if col != col_name]

        df = pd.read_csv(f_in, sep='\t', dtype='str').filter(colmap.keys())
        df.columns = colmap.values()

        # convert to abs floats
        df[col_name] = np.abs(df[col_name].astype(float))

        # Let's filter out duplicates
        duplicates = df.groupby(index_cols)
        df_filtered = pd.DataFrame()

        for _, group in duplicates:
            if len(group) == 1:
                # one pred per variant!
                df_filtered = pd.concat([df_filtered, group], sort=False)
            else:
                # take the biggest absolute value
                df_filtered = pd.concat([df_filtered, group.loc[np.abs(group[col_name]).idxmax()].to_frame().T])

        df_filtered['POS'] = df_filtered['POS'].astype(int)
        df_filtered[col_name] = df_filtered[col_name].astype(float)

        predictors.add(col_name)

        # join with main predictions file
        return predictions.reset_index().merge(df_filtered, on=index_cols, how='left', validate='one_to_one').set_index('ID')

    predictions_mms = helpers.vcf_to_df(os.path.join('variants', 'variants_mmsplice.vcf'))
    predictions_mms = append_from_mms_tsv(predictions_mms, 'pathogenicity', 'MMSplice (Pathogenicity)')
    predictions_mms = append_from_mms_tsv(predictions_mms, 'splicing_efficiency', 'MMSplice (Splicing Efficiency)')

    # join by ID with the other predictions
    df = df.join(predictions_mms.drop(columns=['#CHROM', 'POS', 'REF', 'ALT']))

    df.to_csv(os.path.join('predictions', 'joined', 'predictions.csv'))

    with open(os.path.join('predictions', 'joined', 'predictors.txt'), 'w') as f:
        f.write(','.join(predictors))


def read():
    df = pd.read_csv(os.path.join('predictions', 'joined', 'predictions.csv')).set_index('ID')

    with open(os.path.join('predictions', 'joined', 'predictors.txt'), 'r') as f:
        predictors = f.readline().split(',')

    return df, predictors

if __name__ == '__main__':
    aggregate()