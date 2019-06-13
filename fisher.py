#!/usr/bin/env python

import csv
import numpy as np
import argparse

def parse(handle):
    rows = csv.reader(handle, delimiter=',')
    domains = next(rows)[1:]
    rows = list(rows)
    proteins = [el[0] for el in rows]
    hits = np.array([el[1:] for el in rows])
    print(proteins)
    return domains, proteins, hits

def fisher(N, K, n, x):
    from scipy.special import binom
    res = float(0)
    for i in range(x, N + 1):
        res += (binom(K, i)*binom(N - K, n - i)) / binom(N, n)
    return res

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    # parser.add_argument('-f1', '--f1', help='First csv input filename', required=True)
    # parser.add_argument('-f2', '--f2', help='First csv input filename', required=True)
    args = parser.parse_args()
    NAME1 = 'blast_hmm_out.csv'  # args.f1
    NAME2 = 'nonblast_hmm_out.csv'  # args.f2
    handle1 = open(NAME1, 'r')
    domains1, proteins1, hits1 = parse(handle1)
    handle2 = open(NAME2, 'r')
    domains2, proteins2, hits2 = parse(handle2)
    hits1, hits2 = hits1.astype(np.float), hits2.astype(np.float)
    N = len(proteins1) + len(proteins2)
    print("N=", len(proteins1), "+", len(proteins2))
    K = len(proteins1)
    print("K=", len(proteins1))
    n1, n2 = int(np.sum(hits1)), int(np.sum(hits2))
    n = n1 + n2
    print("n=", n1, "+", n2)
    x = int(np.sum(hits1))
    print("x=", n1)
    res = fisher(N, K, n, x)
    print("Res = ", res)
