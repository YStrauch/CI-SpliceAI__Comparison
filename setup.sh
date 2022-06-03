#!/bin/bash
source src/helpers.sh


if [ -z "$ENV_INSTALLED" ]; then
    conda create --yes -n $ENV python=3.8
fi
conda activate $ENV
conda install --yes numpy pandas matplotlib samtools maxentpy biopython scikit-learn tabulate -c bioconda -c conda-forge
pip install spliceai[cpu] kipoi pyfaidx cispliceai[cpu]

mkdir -p third-party
cd third-party

# -- DOWNLOAD REFERENCE GENOME --
mkdir -p hg
cd hg
curl http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.p13.genome.fa.gz > hg38.fa.gz
gzip -d hg38.fa.gz
cd ..

# -- INSTALL SQUIRLS --
mkdir -p squirls

# get CLI
cd squirls
curl -L https://github.com/TheJacksonLaboratory/Squirls/releases/download/v1.0.0/squirls-cli-1.0.0-distribution.zip > squirls.zip
unzip squirls.zip
rm squirls.zip

# download jannovar
mkdir -p jannovar
cd jannovar
wget ftp://squirls.ielis.xyz/jannovar_v0.35.zip
unzip jannovar_v0.35.zip
rm jannovar_v0.35.zip
cd ..

# download SQUIRLS database
mkdir -p 2103_hg38
cd 2103_hg38
wget ftp://squirls.ielis.xyz/2103_hg38.zip
unzip 2103_hg38.zip
rm 2103_hg38.zip
cd ..

cd ..

# -- INSTALL MMSPLICE THROUGH KIPOI --
mkdir -p mmsplice
cd mmsplice

# download gencode annotations for mmsplice
curl http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz > gencode.gtf.gz
gunzip gencode.gtf.gz
cd ..

# install mmsplice
kipoi_env_name="kipoi-MMSplice__splicingEfficiency"


# create MMSplice environment, if not already exist
installed=$(conda info --envs | grep "$kipoi_env_name " || true;)
if [ -z "$installed" ]; then
    kipoi env create MMSplice/splicingEfficiency
fi
source activate $kipoi_env_name
# install MMSplice into this environment
kipoi env install MMSplice/splicingEfficiency

conda deactivate

cd ..

# pre-process mmsplice gencode
conda activate $kipoi_env_name
python src/setup/mmsplice.py
conda deactivate