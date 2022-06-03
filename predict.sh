#!/bin/bash
source src/helpers.sh

curdir=$(pwd)
HG=$curdir/third-party/hg/hg38.fa

# --- spliceai ---
spliceai -I variants/variants.vcf -O predictions/spliceai.vcf -R $HG -A grch38 -D 4999

# --- ci-spliceai ---
cis-vcf --all -i variants/variants.vcf -o predictions/cis.vcf -a grch38 -d 5000 $HG

# --- ci-spliceai (trained on TRAIN chroms only) ---
python src/predict/cis_train.py


# --- squirls ---
java -jar third-party/squirls/squirls-cli-1.0.0/squirls-cli-1.0.0.jar annotate-vcf src/predict/squirls-config.yml third-party/squirls/jannovar/hg38_refseq.ser variants/variants.vcf predictions/squirls -f vcf

# --- mmsplice ---

# splicing efficiency
conda activate kipoi-MMSplice__splicingEfficiency
cd ~/.kipoi/models/MMSplice/splicingEfficiency
kipoi predict MMSplice/splicingEfficiency -m --dataloader_args='{'gtf': '$curdir/third-party/mmsplice/protein_coding_no_dups.gtf', 'fasta_file': '$HG', 'vcf_file': '$curdir/variants/variants_mmsplice.vcf'}' -o $curdir/predictions/mmsplice_splicing_efficiency.tsv

# pathogenicity
cd ~/.kipoi/models/MMSplice/pathogenicity 
kipoi predict MMSplice/pathogenicity -m --dataloader_args='{'gtf': '$curdir/third-party/mmsplice/protein_coding_no_dups.gtf', 'fasta_file': '$HG', 'vcf_file': '$curdir/variants/variants_mmsplice.vcf'}' -o $curdir/predictions/mmsplice_pathogenicity.tsv
cd $curdir

# --- MES ---
# MES via VEP
cp variants/variants.vcf ~/vep_data/variants_splicing_comparison.vcf
docker run -t -i -v ~/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep vep --distance 5000 --hgvs --vcf --force_overwrite --plugin MaxEntScan,/opt/vep/.vep/Plugins/maxentscan --species homo_sapiens --transcript_version --input_file /opt/vep/.vep/variants_splicing_comparison.vcf  --output_file /opt/vep/.vep/output_splicing_comparison.vcf --database
cp ~/vep_data/output_splicing_comparison.vcf predictions/mes_vep.vcf

# MES sliding window
python src/predict/mes_sliding.py