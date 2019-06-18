#!/usr/bin/env python

import csv
import numpy as np
import argparse
from scipy.special import binom

def parse(handle):
    rows = csv.reader(handle, delimiter=',')
    domains = next(rows)[1:]
    rows = list(rows)
    proteins = [el[0] for el in rows]
    hits = np.array([el[1:] for el in rows])
    print(proteins)
    return domains, proteins, hits

def fisher(N, K, n, x):
    res = float(0)
    for i in range(x, N + 1):
        res += (binom(K, i)/binom(N, n) * binom(N - K, n - i))
    return res

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-f1', '--f1', help='First csv input filename', required=True)
    parser.add_argument('-f2', '--f2', help='Second csv input filename', required=True)
    args = parser.parse_args()
    handle1 = open(args.f1, 'r')
    domains1, proteins1, hits1 = parse(handle1)
    handle2 = open(args.f2, 'r')
    domains2, proteins2, hits2 = parse(handle2)
    hits1, hits2 = hits1.astype(np.float), hits2.astype(np.float)

    N = len(proteins1) + len(proteins2)
    print("N=", len(proteins1), "+", len(proteins2))
    K = len(proteins1)
    print("K=", len(proteins1))

    from scipy.stats import fisher_exact
    # for i, dom in enumerate(domains1):
    #     print('########## - ', dom)
    #     n1 = int(np.sum(hits1[:, i]))
    #     n2 = 0
    #     try:
    #         j = domains2.index(dom)
    #         n2 = int(np.sum(hits2[:, j]))
    #     except ValueError:
    #         pass
    #     n = n1 + n2
    #     # print("n=", n1, "+", n2)
    #     x = n1
    #     # print("x=", x)
    #     res = fisher(N, K, n, x)
    #     # print("Res = ", res)
    #     odds2, res2 = fisher_exact([[x, K - x],[n2, N - K - n2]], 'greater')
    #     print("Res2 = ", res2)

    print("-----------------------")
    results = []
    K = len(proteins2)
    for i, dom in enumerate(domains2):
        print('########## - ', dom)
        n1 = int(np.sum(hits2[:, i]))
        n2 = 0
        try:
            j = domains1.index(dom)
            n2 = int(np.sum(hits1[:, j]))
        except ValueError:
            pass
        n = n1 + n2
        # print("n=", n1, "+", n2)
        x = n1
        # print("x=", x)
        res = fisher(N, K, n, x)
        print("Res = ", res)
        odds2, res2 = fisher_exact([[x, K - x],[n2, N - K - n2]], 'greater')
        print("Res2 = ", res2)
        results.append([dom, str(res), str(res2)])
    reshandle = open('results.csv', 'w')
    a = csv.writer(reshandle, delimiter=',')
    a.writerows(results)

