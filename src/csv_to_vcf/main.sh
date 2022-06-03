#!/bin/bash
source ../helpers.sh

# this script generates the variants.vcf and variants_mmsplice.vcf in the main directory.
# you don't need to run it as the variants.vcf is checked into git.
# if you want to run it, cd into this directory first.

# it parses the variants.csv file into a vcf file. it also normalises (rows and align-left) it for MMSplice, as suggested in the MMSplice pipeline.

HG=../third-party/hg/hg38.fa

# sanity check and conversion to tab-delimited format (plain.vcf)
python convert.py $HG

# Index VCF
bgzip -c plain.vcf > plain.vcf.gz
rm plain.vcf
tabix -fp vcf plain.vcf.gz

# normalise rows, this writes headers as well
bcftools norm -m-both -o norm1.vcf.gz plain.vcf.gz

# move variants into main folder; this is our main file
gzip -d norm1.vcf.gz
cp norm1.vcf ../variants.vcf

# the rest is only needed for MMSplice, see https://github.com/kipoi/models/tree/master/MMSplice
bgzip -c norm1.vcf > norm1.vcf.gz

# normalise left
bcftools norm -f $HG -o norm2.vcf.gz norm1.vcf.gz

# and unzip it again
gzip -d norm2.vcf.gz

# move into main folder
mv norm2.vcf ../variants_mmsplice.vcf

# clean up
rm norm1.vcf
rm norm1.vcf.gz
rm plain.vcf.gz
rm plain.vcf.gz.tbi