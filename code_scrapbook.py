"""
A place to put odd bits of code under construction/for use in notebooks.

Example:
python code_scrapbook.py
"""

import numpy as np
import seaborn as sns
import pybedtools

a = pybedtools.example_bedtool('a.bed')
b = pybedtools.example_bedtool('b.bed')

iap_smush = 

print(a.intersect(b))
