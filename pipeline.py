#!/usr/bin/env python
import os


# os.system('./extend.py')

IN1 = 'input-z2.fasta'
OUT1 = 'nonblast_hmm_out.csv'

# os.system('./scan_pfam.py -i {} -o {}'.format(IN1, 'kurwa'))

IN2 = 'blast_out.fa'
OUT2 = 'blast_hmm_out.csv'
# os.system('./scan_pfam.py -i {} -o {}'.format(IN2, OUT2))

F1 = 'nonblast_hmm_out.csv'
F2 = 'blast_hmm_out.csv'
# os.system('./fisher.py -f1 {} -f2 {}'.format(F1, F2))


import matplotlib.pyplot as plt
import csv

a = csv.reader(open('results.csv'))
assert False, list(a)


