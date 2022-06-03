import os
from cispliceai.model import CISpliceAI
from cispliceai.annotation import Annotator

# load CIS but substitute their models which those from src directory which were trained holding out chromosomes 1,3,5,7,9
cis_train = CISpliceAI(os.path.join('src', 'CI-SpliceAI-TRAIN.pb'))
annotator = Annotator(model=cis_train)

annotator.annotate_vcf(
    reference_path=os.path.join('third-party', 'hg', 'hg38.fa'),
    vcf_in=os.path.join('variants', 'variants.vcf'),
    vcf_out=os.path.join('predictions', 'cis_train.vcf'),
    annotation_table='grch38',
    max_dist_from_var=5000,
    most_significant_only=False
)