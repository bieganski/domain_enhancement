#!/usr/bin/env python
import os


os.system('./extend.py')

IN1 = 'input-z2.fasta'
OUT1 = 'nonblast_hmm_out.csv'

os.system('./scan_pfam.py -i {} -o {}'.format(IN1, 'kurwa'))

IN2 = 'blast_out.fa'
OUT2 = 'blast_hmm_out.csv'

os.system('./scan_pfam.py -i {} -o {}'.format(IN2, OUT2))

F1 = 'nonblast_hmm_out.csv'
F2 = 'blast_hmm_out.csv'

os.system('./fisher.py -f1 {} -f2 {}'.format(F1, F2))


import matplotlib.pyplot as plt
import csv


a = csv.reader(open('results.csv'))
a = list(a)
a = [(name, float(impl1), float(impl2)) for name, impl1, impl2 in a]


# DRAWS PLOT OF DIFFERENCES BETWEEN FISHERS TEST IMPLEMENTATIONS
delta = [abs(impl1 - impl2) for _, impl1, impl2 in a]
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
ax.set_xlabel('Domain')
ax.set_ylabel('Delta')
ax.set_title('Value differences between scipy.stats and custom Fishers test implementation')
plt.plot(delta, 'o')
# plt.show()
plt.savefig('deltas.png')

# DRAWS PLOT OF FISHRS TEST RESULTS FOR EACH DOMAIN
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(111)
ax.set_xlabel('Domain')
ax.set_ylabel('FT')

a = [(dom, ft) for dom, ft, _ in a]
a.sort(key=lambda x: x[1])

x = [dom for dom, _, in a]
y = [ft for _, ft, in a]
plt.plot(x, y, 'o')
plt.xticks(x, x, rotation='vertical')
plt.gcf().subplots_adjust(bottom=0.35)
# plt.show()
plt.savefig('ft.png')

