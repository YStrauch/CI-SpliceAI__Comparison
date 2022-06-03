import cispliceai.fasta
import pysam
import numpy as np
import pandas as pd

fasta = cispliceai.fasta.PyFaidXFasta()

def extract_sequences(genome_path, chromosomes, starts, ends, context=0):
   return [extract_sequence(genome_path, chrom, start, end, context=context) for chrom, start, end in zip(chromosomes, starts, ends)]

def extract_sequence(genome_path, chrom, start, end, context=0):
   start -= context
   end += context
   
   if context:
      # ensure that the reference length is dividable by 2 or delta position will be half R
      end[(end - start)%2 != 0] +=1
     
   return fasta.extract(genome_path, chrom, start, end-start)

def reverse_complement(seq: str):
   return seq.upper().replace('A', 't').replace('C', 'g').replace('G', 'c').replace('T', 'a').upper()[::-1]


def vcf_to_df(file):
   '''Parses a VCF file into a pandas DataFrame. Ignores INFO field.'''
   with pysam.VariantFile(file) as f:
      lines = [line for line in f]
   
   assert np.all([len(line.alts) == 1 for line in lines]), 'More than one ALT per line given'

   df = pd.DataFrame({
      '#CHROM': [line.chrom for line in lines],
      'POS': [line.pos for line in lines],
      'REF': [line.ref for line in lines],
      'ALT': [line.alts[0] for line in lines],
   }, index=[line.id for line in lines])
   df.index.name = 'ID'

   return df