"""
A place to put odd bits of code under construction/for use in notebooks.

Example:
python code_scrapbook.py
"""

import numpy as np
import seaborn as sns
import pybedtools
from pyfaidx import Fasta

a = pybedtools.example_bedtool('a.bed')
b = pybedtools.example_bedtool('b.bed')

iap = pybedtools.BedTool('data/ME_beds/mm10.IAP.mended.meta_subelements.bed')

mm_fasta = pybedtools.BedTool('data/GRCm38.p6.genome.fa')

iap = iap.sequence(fi=mm_fasta)

iap.save_seqs('data/iap_sequence.fa')
